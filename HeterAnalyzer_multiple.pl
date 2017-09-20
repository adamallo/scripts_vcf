#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;
use Env;
#use Parallel::Loops;

##Configuration variables
######################################################

our $sep_param="#";
our $sep_value=":";
our $OFS=",";
our $FS=",";
our $output_vcfs=1;
our $output_list=1;
our $output_comprehensive=1;
our $n_cores=1;
our $covB_sh="covB.sh";
#our $covB_sh="covBdoublecheck.sh";

######################################################

##IO Variables
######################################################
my $execond_inputfile="";
my $filtercond_inputfile="";
my $NABfiltercond_inputfile1="";
my $NABfiltercond_inputfile2="";
my $covBfiltercond_inputfile="";
my $output_file="";
my $output_folder="";
my $Aname="A";
my $Bname="B";
my $Nname="N";
my $SCRIPTSVCF_DIR=$ENV{'SCRIPTSVCF_DIR'};

#Flags
my $help;
my $usage="Usage: $0 [options] -o output_file \n\nWARNING: This is a secondary script that is not inteded to be executed directly.\n\n\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the NAB sample\n\t--NABfilt_cond_inputfile2 : secondary input file for a secondary filter (OR filter implemented in a dirty way) of the NAB sample\n\t--covaltB_cond_inputfile : input file for the filtering taking into account characteristics of the unfiltered in the comparison\n\t--output_dir : output directory for vcf files\n\t--output_folder : output folder\n\t-a/--Aname : Name of the A file\n\t-b/--Bname : Name of the B file\n\t-n/--Nname : Name of the N file\n\t--n_cores : number of cores to execute some steps in parallel (requires the perl package Parallel::Loops)\n\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
######################################################
(! GetOptions(
    'exec_cond_inputfile|e=s' => \$execond_inputfile,
	'filt_cond_inputfile|f=s' => \$filtercond_inputfile,
	'NABfilt_cond_inputfile=s' => \$NABfiltercond_inputfile1,
    'NABfilt_cond_inputfile2=s' => \$NABfiltercond_inputfile2,
    'covaltB_cond_inputfile=s' => \$covBfiltercond_inputfile,
    'output_folder=s' => \$output_folder,
    'Aname|a=s' => \$Aname,
    'Bname|b=s' => \$Bname,
    'Nname|n=s' => \$Nname,
	'output_file|o=s' => \$output_file,
    'n_cores=i' => \$n_cores,
    'help|h' => \$help,
                )) or (($output_file eq "") || ($output_folder eq "") || $help) and die $usage;

##Load Parallel::Loops if it is available and it's needed
#########################################################

if($n_cores>1)
{
    eval "use Parallel::Loops";
    if($@)
    {
        print "\n\nWARNING: You are asking to execute this script using $n_cores cores, but the required module \"Parallel::Loops\" has not been found in \@INC\n\n";
         $n_cores=1;
   }
    else
    {
        print "\nUsing Parallel::Loops with $n_cores cores\n\n";
    }
}

##Input file parsing and directory creation
######################################################

my @exe_parameters=("input");
my @exe_param_values=([("")]);
my @filtering_parameters=("");
my @NABfiltering_parameters1=("");
my @NABfiltering_parameters2=("");
my @covBfiltering_parameters=("");

my @filtering_param_values=([("")]);
my @NABfiltering_param_values1=([("")]);
my @NABfiltering_param_values2=([("")]);
my @covBfiltering_param_values=([("")]);

##Input files

if ($execond_inputfile ne "")
{
    @exe_parameters=();
    parse_parameters_values($execond_inputfile,\@exe_parameters,\@exe_param_values);
}

if ($filtercond_inputfile ne "")
{
    @filtering_parameters=();
    parse_parameters_values($filtercond_inputfile,\@filtering_parameters,\@filtering_param_values);
}

if ($NABfiltercond_inputfile1 ne "")
{
    @NABfiltering_parameters1=();
    parse_parameters_values($NABfiltercond_inputfile1,\@NABfiltering_parameters1,\@NABfiltering_param_values1);
}

if ($NABfiltercond_inputfile2 ne "")
{
    @NABfiltering_parameters2=();
    parse_parameters_values($NABfiltercond_inputfile2,\@NABfiltering_parameters2,\@NABfiltering_param_values2);
}

if ($covBfiltercond_inputfile ne "")
{
    @covBfiltering_parameters=();
    parse_parameters_values($covBfiltercond_inputfile,\@covBfiltering_parameters,\@covBfiltering_param_values);
}

my $vcf_filt_exe="$SCRIPTSVCF_DIR/vcf_filtering.pl";

if (! -f "$vcf_filt_exe")
{
    die "The executable vcf_filtering.pl can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR indicates the directory with this package of scripts\n";
}

## Parsing static variants
############################################################################################

my %N;
my ($ref_n,$nameN);

if ( -f "$Nname.vcf")
{
    ($ref_n,$nameN)=parse_vcf_name("$Nname.vcf");
    #print("DEBUG: N genotype name $nameN\n");
    %N=%{$ref_n};
}
else
{
    die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the variant calling step has finished succesfully and try to execute this script again\n";
}


#Geting NAB
##########################################################################

##Locate NAB
my $pos_N;
my $NABname="";
if( -f "${Nname}_${Aname}_${Bname}.vcf")
{
    $NABname="${Nname}_${Aname}_${Bname}";
    $pos_N=get_col_name_vcf("$NABname.vcf",$nameN);
}
elsif ( -f "${Nname}_${Bname}_${Aname}.vcf")
{
    $NABname="${Nname}_${Bname}_${Aname}";
    $pos_N=get_col_name_vcf("$NABname.vcf",$nameN);
}
else
{
    die "NAB file of $Nname, $Aname and $Bname was not located in the current folder\n";
}

#print("DEBUG: N genotype column $pos_N\n");

@NABfiltering_param_values1=filtervalues_changepost(\@NABfiltering_param_values1,$pos_N);
@NABfiltering_param_values2=filtervalues_changepost(\@NABfiltering_param_values2,$pos_N);

## Main conditions loop
#######################

my $name_condition;
my @filtering_conditions;
our @NABfiltering_conditions1; ##In order to be easily accesed from the filtering functions
our @NABfiltering_conditions2;
our @NABfiltering_conditions;
my @covBfiltering_conditions;
my @exe_conditions;
my %results;
my $bamfiles;
my $parallel;

##Generate all combinations of proposed values for execution and filtering parameters
combs(0,"",\@exe_parameters,\@exe_param_values,\@exe_conditions);
combs(0,"",\@filtering_parameters,\@filtering_param_values,\@filtering_conditions);
combs(0,"",\@NABfiltering_parameters1,\@NABfiltering_param_values1,\@NABfiltering_conditions1);
combs(0,"",\@NABfiltering_parameters2,\@NABfiltering_param_values2,\@NABfiltering_conditions2);
combs(0,"",\@covBfiltering_parameters,\@covBfiltering_param_values,\@covBfiltering_conditions);

#print("DEBUG: @exe_parameters, @exe_param_values");
#print("DEBUG: @exe_conditions,@filtering_conditions,@NABfiltering_conditions1,@NABfiltering_conditions2\n");

if ($n_cores>1)
{
    $parallel = Parallel::Loops->new($n_cores);
    $parallel->share(\%results);
}


## NAB filtering and parsing
##########################################################################################

if($n_cores>1)
{
    $parallel->foreach(\@NABfiltering_conditions1,\&filterNAB);
    $parallel->foreach(\@NABfiltering_conditions2,\&filterNAB);
}
else
{
    foreach my $filtering_condition (@NABfiltering_conditions1)
    {
        filterNAB($filtering_condition);
    }
    foreach my $filtering_condition (@NABfiltering_conditions2)
    {
        filterNAB($filtering_condition);
    }

}
my @NAB_hash_pointers;
my @NAB_hash_pointers1;
my @NAB_hash_pointers2;

for (my $i=0; $i<scalar @NABfiltering_conditions1; ++$i)
{
    $NAB_hash_pointers1[$i]=parse_vcf("$NABname$sep_param$NABfiltering_conditions1[$i].vcf");
}

for (my $i=0; $i<scalar @NABfiltering_conditions2; ++$i)
{
    $NAB_hash_pointers2[$i]=parse_vcf("$NABname$sep_param$NABfiltering_conditions2[$i].vcf");
}

my $pos=0;
for (my $i=0; $i<scalar @NABfiltering_conditions1; ++$i)
{  
    for (my $j=0; $j< scalar @NABfiltering_conditions2; ++$j)
    {
        $NABfiltering_conditions[$pos]="$NABfiltering_conditions1[$i]$sep_param$NABfiltering_conditions2[$j]";
        $NAB_hash_pointers[$pos]=combine_vcf_hashes($NAB_hash_pointers1[$i],$NAB_hash_pointers2[$j]);
        $pos+=1;
    }
}

##Output directory
##################

##Adding original directory before changing to the new one in 
my $dir=getcwd;
my $AfullName="$dir/$Aname";
my $BfullName="$dir/$Bname";
my $NfullName="$dir/$Nname";
my $NABfullName="$dir/$NABname";
mkdir($output_folder);
chdir($output_folder);


##Main Loop
###########
foreach my $exe_condition (@exe_conditions) ##Options that require to call variants again 
{
#print("DEBUG: Exe condition loop $exe_condition\n");
    my @current_conditions=@filtering_conditions;
    my $icond=0;

    if (scalar @covBfiltering_conditions == 0)
    {
        for (my $i=0; $i<scalar @filtering_conditions;++$i)
        {
            $current_conditions[$icond]=[$exe_condition,$filtering_conditions[$i],$covBfiltering_conditions[0]]; ##If we are not using covB this will be still used but undefined
            $icond+=1;
        }
    }
    else
    {
       for (my $i=0; $i<scalar @filtering_conditions;++$i)
       {
           for (my $j=0; $j<scalar @covBfiltering_conditions; ++$j)
           {
               $current_conditions[$icond]=[$exe_condition,$filtering_conditions[$i],$covBfiltering_conditions[$j]];
               $icond+=1;
           }
       }
    }

    #I do not need to paralelize per covBfiltering option, since the different options can be obtained from the same tsv, just considering the variants of interests, if we start with the unfiltered versions (supersets of the others)
    if((! -f "$AfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv") || (! -f "$BfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv"))
    {
       `$covB_sh $AfullName$sep_param$exe_condition.vcf $BfullName$sep_param$exe_condition.vcf $AfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv $BfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv`; 
    }
   
    if($n_cores>1)
    {
        $parallel->foreach(\@current_conditions,\&filter);
    }
    else
    {
        foreach my $filtering_condition (@current_conditions)
        {
            filter($filtering_condition);
        }
    }
}

## Output
#############################################################
open(my $OFILE,">$output_file");
print($OFILE "Sample,Condition,A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean\n");
#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean
my $sample=$output_file;
$sample=basename($sample);
$sample=~s/(.*)\..*/$1/;
foreach my $condition (keys %results)
{
    print($OFILE "$sample$OFS$condition$OFS",array_to_string(@{$results{$condition}}),"\n");
}
close($OFILE);

## MAIN BODY AS FUNCTION FOR PARALLELISM
########################################
sub filter
{
    my $exe_condition;
    my $filtering_condition;
    my $covBfiltering_condition;
    if (defined $_)
    {
        ($exe_condition,$filtering_condition,$covBfiltering_condition)=@{$_};
    }
    else
    {
        ($exe_condition,$filtering_condition,$covBfiltering_condition)=@{$_[0]};
    }
    
    my $condition="$exe_condition$sep_param$filtering_condition";
    my $NcovBcondition="$exe_condition$sep_param$filtering_condition$sep_param$covBfiltering_condition";

    my $filtering_command="$vcf_filt_exe ";
    $filtering_command.=join(" ",split("$sep_value",join(" ",split("$sep_param",$filtering_condition))));
   
    my(%A,%B);
    my @nofilt_results;
    
    if ( -f "$AfullName$sep_param$exe_condition.vcf" && -f "$BfullName$sep_param$exe_condition.vcf")
    {
        %A=%{parse_vcf("$AfullName$sep_param$exe_condition.vcf")};
        %B=%{parse_vcf("$BfullName$sep_param$exe_condition.vcf")};
        my %AN=%{vcf_prune_single(\%A,\%N)};
        my %BN=%{vcf_prune_single(\%B,\%N)};
        @nofilt_results=(scalar keys %A, scalar keys %B, scalar keys %N, scalar keys %AN, scalar keys %BN);
    }
    else
    {
        die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the variant calling step has finished succesfully and try to execute this script again\n";
    }
    
    my %Afilt=%{parse_vcf("$AfullName$sep_param$condition.vcf")};
    my %Bfilt=%{parse_vcf("$BfullName$sep_param$condition.vcf")};

    ##Data to filter variants according to covB. I implemented this in this script instead of in an additional one. From an architectonic point of view may not be ideal, but this will be very quick and don't want to complicate more and more the execution tree of this pipeline.

    my %dataAfiltcovB;
    my %dataBfiltcovB;

    if ( -f "$AfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv" && -f "$BfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv")
    {
        %dataAfiltcovB=%{parse_tsv("$AfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv")};
        %dataBfiltcovB=%{parse_tsv("$BfullName$sep_param$exe_condition${sep_param}covBfiltering.tsv")};
    }
    else
    {
        die "Missing tsv files. There may be an error in the internal reference to \$covB_sh\n";
    }

# Right now I'm keeping a lot of hashes in memory instead of reusing variables. I do it just in case I need them in posterior statistics/calculations etc.
# We could be interested on changing this if the performance were really bad.

#Compare Afilt with B without filter --> Common variants + %
    my @statsAfilt;
    my ($ref_common_variantsAfilt,$ref_different_variantsAfilt)=vcf_compare_parsed(\%B,\%Afilt,\@statsAfilt); ##I have to generate two hashes. One with common variants, the other with non common. Thus, the consecutive filter I can do it towards these new (smallest) hashes.
#Compare Bfiltered with A without filter --> Common variants + %
    my @statsBfilt;
    my ($ref_common_variantsBfilt,$ref_different_variantsBfilt)=vcf_compare_parsed(\%A,\%Bfilt,\@statsBfilt);

#Stats filter
    my @statsfiltU;
    my @statsfiltI;
    my ($ref_common_variantsfiltU,$ref_different_variantsfiltU)=vcf_unite_parsed($ref_common_variantsAfilt,$ref_different_variantsAfilt,$ref_common_variantsBfilt,$ref_different_variantsBfilt,\@statsfiltU);
    my ($ref_common_variantsfiltI,$ref_different_variantsfiltI)=vcf_intersect_parsed($ref_common_variantsAfilt,$ref_different_variantsAfilt,$ref_common_variantsBfilt,$ref_different_variantsBfilt,\@statsfiltI);
    my @statsfiltmean=(($statsAfilt[0]+$statsBfilt[0])/2.0,($statsAfilt[1]+$statsBfilt[1])/2.0,($statsAfilt[2]+$statsBfilt[2])/2.0);

#Substract N from Afilt. Compare the result to B without filter --> Common variants + %
    my @statsAfiltN;
    my ($ref_common_variantsAfiltN,$ref_different_variantsAfiltN)=vcf_prune($ref_common_variantsAfilt,$ref_different_variantsAfilt,\%N,\@statsAfiltN);

#Substract N from Bfilt. Compare the result to A without filter --> Common variants + %
    my @statsBfiltN;
    my ($ref_common_variantsBfiltN,$ref_different_variantsBfiltN)=vcf_prune($ref_common_variantsBfilt,$ref_different_variantsBfilt,\%N,\@statsBfiltN);

#Mean stats filterN
    my @statsfiltNU;
    my @statsfiltNI;
    my ($ref_common_variantsfiltNU,$ref_different_variantsfiltNU)=vcf_prune($ref_common_variantsfiltU,$ref_different_variantsfiltU,\%N,\@statsfiltNU);
    my ($ref_common_variantsfiltNI,$ref_different_variantsfiltNI)=vcf_prune($ref_common_variantsfiltI,$ref_different_variantsfiltI,\%N,\@statsfiltNI);

    my @statsfiltNmean=(($statsAfiltN[0]+$statsBfiltN[0])/2.0,($statsAfiltN[1]+$statsBfiltN[1])/2.0,($statsAfiltN[2]+$statsBfiltN[2])/2.0);

#CovB
    
    my @statsAfiltNcovB;
    my @statsBfiltNcovB;
    my @statsfiltNcovBU;
    my @statsfiltNcovBI;
    
    #AfiltNcovB
    my ($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB)=vcf_prune_covB($ref_common_variantsAfiltN,$ref_different_variantsAfiltN,\%dataBfiltcovB,\@statsAfiltNcovB); 

    #BfiltNcovB
    my ($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB)=vcf_prune_covB($ref_common_variantsBfiltN,$ref_different_variantsBfiltN,\%dataAfiltcovB,\@statsBfiltNcovB);

    #Union 
    my ($ref_common_variantsfiltNcovBU,$ref_different_variantsfiltNcovBU)=vcf_unite_parsed($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,\@statsfiltNcovBU);

    #Intersection
    my ($ref_common_variantsfiltNcovBI,$ref_different_variantsfiltNcovBI)=vcf_intersect_parsed($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,\@statsfiltNcovBI);
    
    #Means
    my @statsfiltNcovBmean=(($statsAfiltNcovB[0]+$statsBfiltNcovB[0])/2.0,($statsAfiltNcovB[1]+$statsBfiltNcovB[1])/2.0,($statsAfiltNcovB[2]+$statsBfiltNcovB[2])/2.0);

#Output of list of variants and/or intermediate vcf files
    if($output_list)
    {
        if($output_comprehensive)
        {
            write_variant_list($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.list");
            write_variant_list($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.list");
            write_variant_list($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.list");
            write_variant_list($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.list");
            write_variant_list($ref_common_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_common.list");
            write_variant_list($ref_common_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_common.list");
            write_variant_list($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_different.list");
            write_variant_list($ref_different_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_different.list");
            write_variant_list($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_different.list");
            write_variant_list($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_different.list");
            write_variant_list($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.list");

        }

        write_variant_list($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.list");
        write_variant_list($ref_common_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_common.list");
        write_variant_list($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.list");
        write_variant_list($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.list");
        write_variant_list($ref_common_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_common.list");
        write_variant_list($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.list");

    }

    if($output_vcfs)
    {
        if($output_comprehensive)
        {
            write_variant_vcf($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_common.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_common.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_different.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_different.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");

        }
        write_variant_2vcf($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");

    } 

    for (my $i=0; $i< scalar @NABfiltering_conditions; ++$i)
    {

        my $NAB_condition=$NABfiltering_conditions[$i];
        my $NAB=$NAB_hash_pointers[$i];

#Substract NAB from AfiltN. Compare the results to B without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
        my @statsAfiltNAB;
        my ($ref_common_variantsAfiltNAB,$ref_different_variantsAfiltNAB)=vcf_prune($ref_common_variantsAfiltN,$ref_different_variantsAfiltN,$NAB,\@statsAfiltNAB);

#Substract NAB from BfiltN. Compare the results to A without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
        my @statsBfiltNAB;
        my ($ref_common_variantsBfiltNAB,$ref_different_variantsBfiltNAB)=vcf_prune($ref_common_variantsBfiltN,$ref_different_variantsBfiltN,$NAB,\@statsBfiltNAB);

#Substract NAB from AfiltNcovB. Compare the results to B without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
        my @statsAfiltcovBNAB;
        my ($ref_common_variantsAfiltcovBNAB,$ref_different_variantsAfiltcovBNAB)=vcf_prune($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$NAB,\@statsAfiltcovBNAB);

#Substract NAB from BfiltNcovB. Compare the results to A without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
        my @statsBfiltcovBNAB;
        my ($ref_common_variantsBfiltcovBNAB,$ref_different_variantsBfiltcovBNAB)=vcf_prune($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,$NAB,\@statsBfiltcovBNAB);

#Mean stats filter NAB		
        my @statsfiltNABU;
        my @statsfiltNABI;

        my ($ref_common_variantsfiltNABU,$ref_different_variantsfiltNABU)=vcf_prune($ref_common_variantsfiltNU,$ref_different_variantsfiltNU,$NAB,\@statsfiltNABU);
        my ($ref_common_variantsfiltNABI,$ref_different_variantsfiltNABI)=vcf_prune($ref_common_variantsfiltNI,$ref_different_variantsfiltNI,$NAB,\@statsfiltNABI);

        my @statsfiltNABmean=(($statsAfiltNAB[0]+$statsBfiltNAB[0])/2.0,($statsAfiltNAB[1]+$statsBfiltNAB[1])/2.0,($statsAfiltNAB[2]+$statsBfiltNAB[2])/2.0);

#Mean stats filter covBNAB		
        my @statsfiltcovBNABU;
        my @statsfiltcovBNABI;

        my ($ref_common_variantsfiltcovBNABU,$ref_different_variantsfiltcovBNABU)=vcf_prune($ref_common_variantsfiltNU,$ref_different_variantsfiltNU,$NAB,\@statsfiltcovBNABU);
        my ($ref_common_variantsfiltcovBNABI,$ref_different_variantsfiltcovBNABI)=vcf_prune($ref_common_variantsfiltNI,$ref_different_variantsfiltNI,$NAB,\@statsfiltcovBNABI);

        my @statsfiltcovBNABmean=(($statsAfiltcovBNAB[0]+$statsBfiltcovBNAB[0])/2.0,($statsAfiltcovBNAB[1]+$statsBfiltcovBNAB[1])/2.0,($statsAfiltcovBNAB[2]+$statsBfiltcovBNAB[2])/2.0);

#Output of list of variants and/or intermediate vcf files
        if($output_list)
        {
            if($output_comprehensive)
            {
                write_variant_list($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
                write_variant_list($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
                write_variant_list($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_common_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
                write_variant_list($ref_common_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
                write_variant_list($ref_different_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list");

            }

            write_variant_list($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
            write_variant_list($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
            write_variant_list($ref_common_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
            write_variant_list($ref_common_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
        }
        if($output_vcfs)
        {
            if($output_comprehensive)
            {
                write_variant_vcf($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_2vcf($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_2vcf($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");

                write_variant_vcf($ref_common_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_common_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_different_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$AfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_different_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$BfullName$sep_param$condition.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_2vcf($ref_different_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_2vcf($ref_different_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition");

            }
            write_variant_2vcf($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
            write_variant_2vcf($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
            write_variant_2vcf($ref_common_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
            write_variant_2vcf($ref_common_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$AfullName$sep_param$condition.vcf","$BfullName$sep_param$condition.vcf","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
        }

#Store and/or print
        my @statistics=(@nofilt_results,@statsAfilt,@statsBfilt,@statsfiltU,@statsfiltI,@statsfiltmean,@statsAfiltN,@statsBfiltN,@statsfiltNU,@statsfiltNI,@statsfiltNmean,@statsAfiltNcovB,@statsBfiltNcovB,@statsfiltNcovBU,@statsfiltNcovBI,@statsAfiltNAB,@statsBfiltNAB,@statsfiltNABU,@statsfiltNABI,@statsfiltNABmean,@statsAfiltcovBNAB,@statsBfiltcovBNAB,@statsfiltcovBNABU,@statsfiltcovBNABI);
#Condition,#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean
        $results{"$NcovBcondition${sep_param}NAB$sep_param$NAB_condition"}=\@statistics;
#print("DEBUG:$NcovBcondition$OFS",array_to_string(@statistics),"\n");
    }
}


###################################################################################
###FUNCTIONS
###################################################################################

#Recursive function to generate parameter combinations
#One recursion per parameter
#########################################################
sub combs 
{
    my ($id,$string,$ref_array1,$ref_array2,$ref_array_output)=@_;
    my $priv_sep_value=$sep_value;
#print("DEBUG: $id,$string,$ref_array1,$ref_array2,$ref_array_output\n");
    my @values=@{@{$ref_array2}[$id]};
    my $parameter=@{$ref_array1}[$id];
    if ($parameter=~/=$/ )
    {
        $priv_sep_value="";
    }   

    if($id<scalar(@{$ref_array1})-1) #If there are more params recurse
    {
        foreach my $value (@values)
        {
            combs($id+1,"$string$parameter$priv_sep_value${value}$sep_param",$ref_array1,$ref_array2,$ref_array_output);
        }
    }
    else #Store the string in the output array
    {
        foreach my $value (@values)
        {
            my $new_string="$string$parameter$priv_sep_value$value";
            push(@{$ref_array_output},$new_string);
        }
    }
}

#Parses parameters file, into an array of parameters
#and a matrix of values per parameter
#####################################################
sub parse_parameters_values
{
    my ($file,$ref_params,$ref_values)=@_;
    open (my $FILE, $file) or die "Error opening the file $file";
    my @content=<$FILE>;
    my $i=0;
    foreach my $line (@content)
    {
        chomp($line);
        my @temp=split($FS,$line);
        push(@{$ref_params},splice(@temp,0,1));
        ${$ref_values}[$i]=\@temp;
        ++$i;
    }
    close($FILE);
}

#Generates an string concatenating array members with OFS
#####################################################
sub array_to_string
{
    my $outstring="";
    for (my $i=0;$i<scalar(@_);++$i)
    {
        $outstring.="$_[$i]$OFS";
    }
    chop($outstring);
    return $outstring;
}

#Parse vcf
###########################################
sub parse_vcf
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            $key=$vcf1[$i];
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/; ####TODO: This may be quicker with a split+join strategy. To check it!
                chomp($key);
            $hash{$key}=1;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return \%hash;

}

#Parse tsv
###########################################
sub parse_tsv
{
    my ($tsv1_file)=@_;
    open(my $TSV1,$tsv1_file) or die "The file $tsv1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @tsv1=<$TSV1>;
    close($TSV1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    my @values;

    for ($i=0;$i<scalar @tsv1;$i++)
    {
        chomp($tsv1[$i]);
        @values=split("\t",$tsv1[$i]);
        splice(@values,0,2);
        $key="$values[0]$OFS$values[1]";
        $hash{$key}=\@values;
        #print("DEBUG: New tsv variant data being hashed $key\n");
    }
    return \%hash;

}

#Parse vcf and get name of the sample (for only one sample)
###########################################################
sub parse_vcf_name
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    my $name;

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $name=(split("\t",$vcf1[$i-1]))[9];
                chomp($name);
#print("DEBUG: Name $name\n");
                $flag=1;
            }
            $key=$vcf1[$i];
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/; ####TODO: This may be quicker with a split+join strategy. To check it!
                chomp($key);
            $hash{$key}=1;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return (\%hash,$name);

}
#Compare two vcf files and generate statistics
##############################################
sub vcf_compare2
{
    my ($vcf1_file,$vcf2_file,$ref_statistics)=@_;

    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    open(my $VCF2,$vcf2_file) or die "The file $vcf2_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf1=<$VCF1>;
    my @vcf2=<$VCF2>;
    close($VCF1);
    close($VCF2);
    clean_vcfcontent(\@vcf1);
    clean_vcfcontent(\@vcf2);
    my %variants1=%{variants_to_hash(\@vcf1)};
    my %variants2=%{variants_to_hash(\@vcf2)};
    my %commonvariants;
    my $n1=scalar(keys %variants1);
    my $n2=scalar(keys %variants2);

    if ($n1<$n2) ##This could be improved in terms on speed following the strategy of the next subroutine. I'm not using this one any more, so I won't update it... at least so far.
    {
        foreach my $variant (keys %variants1)
        {
            if (exists $variants2{$variant})
            {
                $commonvariants{$variant}=1;
#print("DEBUG: New common variant $variant");
            }
        }
    }
    else
    {
        foreach my $variant (keys %variants2)
        {
            if (exists $variants1{$variant})
            {
                $commonvariants{$variant}=1;
#print("DEBUG: New common variant $variant");
            }
        }
    }
    @{ $ref_statistics }=($n1,$n2,scalar (keys %commonvariants));
    return \%commonvariants;

}

#Compare two variant hashes, return a hash with commonvariants and another one with variants in the problem
# that aren't present in the reference and generate statistics
##############################################
sub vcf_compare_parsed
{
    my ($vcf_reference,$vcf_problem,$ref_statistics)=@_;
    my %variants2=%{$vcf_problem};
    my %commonvariants;
    my $n1=scalar(keys %{$vcf_reference});
    my $n2=scalar(keys %variants2);


    foreach my $variant (keys %{$vcf_reference})
    {
        if (exists $variants2{$variant})
        {
            $commonvariants{$variant}=1;
            delete $variants2{$variant}; ###This is a local copy of the original hash. I cut it down to reduce the number of comparisons.
#print("DEBUG: New common variant $variant");
        }

        if(scalar(keys %variants2)== 0)
        {
            #print("DEBUG: Variants2 is empty\n");
            last; ##No more pending comparisons
        }
    }	

    if ($n1<$n2)
    {
        print("WARNING: The number of variants in the unfiltered sample is smaller than the one of the filtered\n");		
    }

    my $n_common=scalar(keys %commonvariants);
    if($n2==0)
    {
        @{ $ref_statistics }=(0,$n_common,$n2); ##Stats= proportion of selected reads in the reference, number of selected variants

    }
    else
    {
        @{ $ref_statistics }=($n_common/$n2,$n_common,$n2); ##Stats= proportion of selected reads in the reference, number of selected variants
    }
    return (\%commonvariants,\%variants2); ##Common, not_in_ref

}

#Filter out variants in one hash from two others and generate statistics
########################################################################
sub vcf_prune
{
    my ($vcf_1,$vcf_2,$vcf_todelete,$ref_statistics)=@_;
    my %variants1=%{$vcf_1};
    my %variants2=%{$vcf_2};
    my $tag1=0;
    my $tag2=0;

    foreach my $variant_to_remove (keys %{$vcf_todelete})
    {
#print("DEBUG: Variant $variant_to_remove ");
		if ($tag1==0 and exists($variants1{$variant_to_remove}))
		{
			delete($variants1{$variant_to_remove});
            #print("deleted in vcf1");
		}
		if ($tag2==0 and exists($variants2{$variant_to_remove}))
		{
			delete($variants2{$variant_to_remove});
            #print("deleted in vcf2");
		}
		#print("\n");
		if($tag1==0 and scalar(keys %variants1)==0)
		{
			$tag1=1;
            #print("DEBUG: No more variants in vcf1\n");
		}	
		if($tag2==0 and scalar(keys %variants2)==0)
		{
			$tag2=1;
            #print("DEBUG: No more variants in vcf2\n");
		}
		
		if($tag1==1 and $tag2==1)
		{
            #print("DEBUG: No more variants\n");
			last;
		}
		
	}
	my $n_selvariants=scalar keys %variants1;
	my $n_filtvariants=scalar keys %variants2;
    if($n_selvariants+$n_filtvariants==0)
    {
         @{ $ref_statistics }=(0,$n_selvariants,$n_selvariants+$n_filtvariants); 
    }
    else
    {   	
	    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
	}
    return (\%variants1,\%variants2);
}

#Filter out variants in one hash depending on certain filters from two others and generate statistics
#####################################################################################################
sub vcf_prune_covB
{
    my ($vcf_1,$vcf_2,$tsv_todelete,$ref_statistics,$covBfiltering_options)=@_;
    my %variants1=%{$vcf_1};
    my %variants2=%{$vcf_2};
    my $tag1=0;
    my $tag2=0;
    my $totalreads=0;
    my $totalaltreads=0;
    my $prop_alts=0;

    #Parser for covBfiltering_options 
    #################################################################################
    
    my $min_coverage=0;
    my $max_alternative=9**9**9; ##pseudo +inf
    my $max_propalt=1;
    
    my @covBfiltering_params=split($sep_param,$covBfiltering_options);
    foreach my $covBoption (@covBfilteriong_params)
    {
        my ($param, $value)= split($sep_value,$covBoption);

        if ($param =~ /--min_coverage/i)
        {
           $min_coverage=$value; 
        }
        elsif ($param =~ /--max_alternative/i)
        {
           $max_alternative=$value;
        }
        elsif ($param =~ /--max_propalt/i)
        {
           $max_propalt=$value;
        }
        else
        {
            die "The param $param has not been recognized properly in the covB filtering step\n";
        }
    }
    print("DEBUG:Final filtering options, minimum coverage in B: $min_coverage, maximum number of reads supporting the alternative allele in B: $max_alternative, maximum proportion of alternative alleles in B: $max_propalt\n");
   
    ##Filtering loop
    ############### 
    my $remove_var=0;

    foreach my $tsv_variant (keys %{$tsv_todelete})
    {

        #print("DEBUG: Variant $variant_to_remove ");
        $totalreads=0;
        $totalaltreads+=$_ for split(",",${$tsv_todelete{$tsv_variant}}[3]); ##Oneline for loop
        $totalreads=$totalaltreads+${$tsv_todelete{$tsv_variant}}[2];
        $prop_alts=($totalaltreads+0.0)/$totalreads;
        
        if($totalreads < $min_coverage || $totalaltreads>$max_alternative || $prop_alts>$max_propalt)
        {
            $remove_var=1;
        } 
        else
        {
            $remove_var=0;
        }
        
        print("DEBUG: totalreads $totalreads, totalalternative $totalaltreads, proportion alternative $prop_alts,var 1 ${$tsv_todelete{$tsv_variant}}[2], var 2 ${$tsv_todelete{$tsv_variant}}[3]\n");###$tsv_todelete[2] Number of bases in the reference allele $tsv_todelete[3]

		if ($tag1==0 and exists($variants1{$variant_to_remove}) and $remove_var==1) 
		{
			delete($variants1{$variant_to_remove});
            print("DEBUG: deleted in vcf1\n");
		}
		if ($tag2==0 and exists($variants2{$variant_to_remove}) and $remove_var==1)
		{
			delete($variants2{$variant_to_remove});
            print("DEBUG: deleted in vcf2\n");
		}
		#print("\n");
		if($tag1==0 and scalar(keys %variants1)==0)
		{
			$tag1=1;
            print("DEBUG: No more variants in vcf1\n");
		}	
		if($tag2==0 and scalar(keys %variants2)==0)
		{
			$tag2=1;
            print("DEBUG: No more variants in vcf2\n");
		}
		
		if($tag1==1 and $tag2==1)
		{
            print("DEBUG: No more variants\n");
			last;
		}
		
	}
	my $n_selvariants=scalar keys %variants1;
	my $n_filtvariants=scalar keys %variants2;
    if($n_selvariants+$n_filtvariants==0)
    {
         @{ $ref_statistics }=(0,$n_selvariants,$n_selvariants+$n_filtvariants); 
    }
    else
    {   	
	    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
	}
    return (\%variants1,\%variants2);
}

#Filter out variants in one hash from another
#############################################
sub vcf_prune_single
{
    my ($vcf_1,$vcf_todelete)=@_;
    my %variants1=%{$vcf_1};

    foreach my $variant_to_remove (keys %{$vcf_todelete})
    {
#print("DEBUG: Variant $variant_to_remove ");
		if (exists($variants1{$variant_to_remove}))
		{
			delete($variants1{$variant_to_remove});
            #print("deleted in vcf1");
		}
		#print("\n");
		if(scalar(keys %variants1)==0)
		{
            last;
            #print("DEBUG: No more variants in vcf1\n");
		}
		
	}
    return (\%variants1);
}

#Combines the variants of two samples. Calculates the union of the two samples
#(union of the commons and union of the differences) and remove spurious
#differences (that can potentially appear in the union).Then calculate the stats.
########################################################################
sub vcf_unite_parsed
{
    my ($ref_common_variantsA,$ref_different_variantsA,$ref_common_variantsB,$ref_different_variantsB,$ref_statistics)=@_;
    my %common_variants=(%{$ref_common_variantsA},%{$ref_common_variantsB}); ##Union
    my %different_variants=(%{$ref_different_variantsA},%{$ref_different_variantsB});#Union

    foreach my $common_variant (keys %common_variants) ##We have to delete possible variants that are not different any more
    {
        delete($different_variants{$common_variant}); 
    }
    
    my $n_selvariants=scalar keys %common_variants;
    my $n_filtvariants=scalar keys %different_variants;
    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    return (\%common_variants,\%different_variants);
}

#Combines the variants of two samples. Calculates the intersection of the two samples
#(intersection of the commons and recalculates the differences from all the variants) and remove spurious
#Then calculate the stats.
########################################################################
sub vcf_intersect_parsed
{
    my ($ref_common_variantsA,$ref_different_variantsA,$ref_common_variantsB,$ref_different_variantsB,$ref_statistics)=@_;
    my %common_variants;
    
    foreach my $key (keys %$ref_common_variantsA)
    {
        if (exists $ref_common_variantsB->{$key})
        {
            $common_variants{$key}=$ref_common_variantsA->{$key};
        }
    }
    
    my %different_variants=(%{$ref_different_variantsA},%{$ref_different_variantsB},%{$ref_common_variantsA},%{$ref_common_variantsB});#Union of all variants

    foreach my $common_variant (keys %common_variants) ##We recalculate the different variants (total-common)
    {
        delete($different_variants{$common_variant});
    }
   
    my $n_selvariants=scalar keys %common_variants;
    my $n_filtvariants=scalar keys %different_variants;
    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    return (\%common_variants,\%different_variants);
}

sub combine_vcf_hashes
{
    my($ref1,$ref2)=@_;
    my %common_variants=(%{$ref1},%{$ref2}); ##Union
    return \%common_variants;
}

#Removes the VCF header
##############################################
sub clean_vcfcontent
{
	my ($array)=@_;
	my @array;
	my $last_comment=0;
	for (my $i=0;$i<scalar @{$array};$i++)
	{
		if(${$array}[$i]=~/^##/)
		{
			++$last_comment;
		}
		else
		{
			last;
		}
	}
	splice(@{ $array },0,$last_comment);
}

# Generates a hash of variants from an array with the content of a VCF without the header
#########################################################################################
sub variants_to_hash
{
	my ($array)=@_;
	my %hash;
	my $key;
	for (my $i=1; $i<scalar @{$array};$i++) ###The first line is the header
	{
		$key=${$array}[$i];
		$key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
		$hash{$key}=1;
		#print("DEBUG: New variant being hashed $key\n");
	}
	return \%hash;
}


# Writes a list of variants in csv format contained in a hash
# ############################################################
sub write_variant_list
{
    my ($ref_hash,$filename)=@_;
    open(my $FILE, ">$filename");
    print($FILE "#CHROM,POS\n");
    foreach my $variant (keys %{$ref_hash})
    {
        print($FILE join($OFS,split(",",$variant)),"\n");
    }
    close $FILE;
}

# Writes a vcf with the variables contained in a hash selected from another VCF file
# ##################################################################################

sub write_variant_vcf
{
    my ($ref_hash,$filename,$vcf,$comment)=@_;
    open(my $OFILE, ">$filename") or die "Error opening the file $vcf\n";
    open(my $IFILE, "$vcf") or die "Error opening the file $vcf\n";
    my @icontent= <$IFILE>;
    close($IFILE);
    my $flag=0;
    my %hash=%{$ref_hash};
    my $key;
     
    #Copying the header and adding a new line with filtering info
    #Then adding the variants that are present in the hash.
    for(my$i=0;$i< scalar @icontent; ++$i)
    {
        if ($flag==0 and $icontent[$i]=~/^##/)
        {
            print($OFILE $icontent[$i]);
        }
        elsif($flag==0)
        {
            print($OFILE "$comment\n$icontent[$i]");
            $flag=1;
        }
        else
        {
            $key=$icontent[$i];
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
            chomp($key);
            #print("DEBUG: Key $key\n");
            if(exists $hash{$key})
            {
                print($OFILE $icontent[$i]);
                delete($hash{$key});
            }
            if(scalar keys %hash == 0)
            {
                last;
            }
            
        }

    }

    close $OFILE;
}

# Writes a vcf with the variables contained in a hash selected from two VCF files
# ATTENTION: Variants in the first VCF file have higher priority when they are 
# present in the two VCF files. The header comes also from the first VCF file.
# ##################################################################################

sub write_variant_2vcf
{
    my ($ref_hash,$filename,$vcf,$vcf2,$comment)=@_;
    open(my $OFILE, ">$filename");
    open(my $IFILE, $vcf);
    open(my $IFILE2, $vcf2);
    my @icontent= <$IFILE>;
    my @icontent2= <$IFILE2>;
    close($IFILE);
    close($IFILE2);
    my $flag=0;
    my %hash=%{$ref_hash};
    my $key;
     
    #Copying the header and adding a new line with filtering info
    #Then adding the variants that are present in the hash.
    for(my $i=0;$i< scalar @icontent; ++$i)
    {
        if ($flag==0 and $icontent[$i]=~/^##/)
        {
            print($OFILE $icontent[$i]);
        }
        elsif($flag==0)
        {
            print($OFILE "$comment\n$icontent[$i]");
            $flag=1;
        }
        else
        {
            $key=$icontent[$i];
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
            chomp($key);
            #print("DEBUG: Key $key\n");
            if(exists $hash{$key})
            {
                print($OFILE $icontent[$i]);
                delete($hash{$key});
            }
            if(scalar keys %hash == 0)
            {
                last;
            }
            
        }

    }
    if (scalar keys %hash !=0)
    {
            for(my $i=0;$i< scalar @icontent2; ++$i)
    	    {
		unless ($icontent2[$i]=~/^#/)
		{
			$key=$icontent2[$i];
            		$key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
            		chomp($key);
            		#print("DEBUG: Key $key\n");
            		if(exists $hash{$key})
            		{
                		print($OFILE $icontent2[$i]);
                		delete($hash{$key});
            		}
            		if(scalar keys %hash == 0)
            		{
                		last;
            		}
		}
	    }
    }

    close $OFILE;
}

## Generates the vcf files for NAB
########################################
sub filterNAB
{
        #my $exe_condition;
        my $filtering_condition;
        if (defined $_)
        {
            #($exe_condition,$filtering_condition)=@{$_};
            $filtering_condition=$_;
        }
        else
        {
            ##($exe_condition,$filtering_condition)=@{$_[0]};
            $filtering_condition=$_[0];
        }
    #my $condition="$exe_condition$sep_param$filtering_condition";
    my $condition="$filtering_condition";
    my $filtering_command="$vcf_filt_exe ";
    $filtering_command.=join(" ",split("$sep_value",join(" ",split("$sep_param",$filtering_condition))));

    if (!-f "$NABname$sep_param$condition.vcf")
    {
        #print("Filtering $NABname$sep_param$exe_condition.vcf to generate $NABname$sep_param$condition.vcf\n");
        print("Filtering $NABname.vcf to generate $NABname$sep_param$condition.vcf\n");
        ## Filter the right vcf file generated in the outside loop
        #system("$filtering_command -i $NABname$sep_param$exe_condition.vcf -o $NABname$sep_param$condition.vcf");
        system("$filtering_command -i $NABname.vcf -o $NABname$sep_param$condition.vcf");
    }
    else
    {
        print("$NABname$sep_param$condition.vcf has been previously generated and it will be recycled\n");
    }
}

## Gets the position in a vcf file of the sample with the given name
####################################################################

sub get_col_name_vcf
{
    my ($vcf,$name)=@_;
    open(my $VCF,$vcf) or die "The file $vcf is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf=<$VCF>;
    close($VCF);
    my $flag=0;
    my $i;
    my $j;
    my @cols;
    my $rescol=-1;;

    for ($i=0;$i<scalar @vcf;$i++)
    {
        unless($flag==0 and $vcf[$i]=~/^#/)
        {
            if ($flag==0)
            {
                        @cols=split("\t",$vcf[$i-1]);
                        for ($j=0;$j<scalar @vcf;$j++)
                        {
                            if($cols[$j]=~/$name/)
                            {
                                $rescol=$j-9; ### VCF files have 9 columns before the first genotype
                                last;
                            }
                        }
                        #print("DEBUG: Name $name\n");
                        $flag=1;
                        last;
            }
        }
    } 
    if ($rescol < 0)
    {
        die "Column with genotype $name not detected \n";
    }
    return $rescol;
}

# Change the genotype column from 0 to the one indicated by $pos. 
# Pos is obtained with parse_vcf_name + get_col_name_vcf
################################################################

sub filtervalues_changepost
{
    my ($ref_array,$pos)=@_;
    my @array=@{$ref_array};
    for (my $i=0; $i< scalar @array; $i++)
    {
        my @temp=@{$array[$i]};
    #    print("DEBUG: before modification ",join(",",@temp));
        for (my $j=0; $j<scalar @temp; $j++)
        {
            $temp[$j]=~s/0_/${pos}_/;
        }
        $array[$i]=\@temp;
    #    print(" after modification ",join(",",@{$array[$i]}),"\n");
    }
    return @array;
   
}
