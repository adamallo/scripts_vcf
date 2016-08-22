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
our $output_list=0;
our $output_comprehensive=1;
our $n_cores=1;

######################################################

##IO Variables
######################################################
my $execond_inputfile="";
my $filtercond_inputfile="";
my $NABfiltercond_inputfile1="";
my $NABfiltercond_inputfile2="";
my $output_file="";
my $SCRIPTSVCF_DIR=$ENV{'SCRIPTSVCF_DIR'};

#Flags
my $help;
my $usage="Usage: $0 [options] -o output_file \n\nWARNING: This is a secondary script that is not inteded to be executed directly.\n\n\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the NAB sample\n\t--NABfilt_cond_inputfile2 : secondary input file for a secondary filter (OR filter implemented in a dirty way) of the NAB sample\n\t--output_dir : output directory for vcf files\n\t--n_cores : number of cores to execute some steps in parallel (requires the perl package Parallel::Loops)\n\n";
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
	'output_file|o=s' => \$output_file,
    'n_cores=i' => \$n_cores,
    'help|h' => \$help,
                )) or (($output_file eq "") || $help) and die $usage;

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
my @filtering_parameters=("");
my @NABfiltering_parameters1=("");
my @NABfiltering_parameters2=("");
my @exe_param_values=([("")]);
my @filtering_param_values=([("")]);
my @NABfiltering_param_values1=([("")]);
my @NABfiltering_param_values2=([("")]);

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

my $vcf_filt_exe="$SCRIPTSVCF_DIR/vcf_filtering.pl";

if (! -f "$vcf_filt_exe")
{
    die "The executable vcf_filtering.pl can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR indicates the directory with this package of scripts\n";
}

## Parsing static variants
############################################################################################

my (%A,%B,%N);
my ($ref_n,$nameN);
my @nofilt_results;

if ( -f "A.vcf" && -f "B.vcf" && -f "N.vcf")
{
    %A=%{parse_vcf("A.vcf")};
    %B=%{parse_vcf("B.vcf")};
    ($ref_n,$nameN)=parse_vcf_name("N.vcf");
    #print("DEBUG: N genotype name $nameN\n");
    %N=%{$ref_n};
    my %AN=%{vcf_prune_single(\%A,\%N)};
    my %BN=%{vcf_prune_single(\%B,\%N)};
    @nofilt_results=(scalar keys %A, scalar keys %B, scalar keys %N, scalar keys %AN, scalar keys %BN);
}
else
{
    die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the variant calling step has finished succesfully and try to execute this script again\n";
}


## Fixing NAB filters, since the N column is not always 0 as expected (weird platypus behaviour)
##########################################################################

my $pos_N=get_col_name_vcf("NAB.vcf",$nameN);

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
my @exe_conditions;
my %results;
my $bamfiles;
my $parallel;

##Generate all combinations of proposed values for execution and filtering parameters
combs(0,"",\@exe_parameters,\@exe_param_values,\@exe_conditions);
combs(0,"",\@filtering_parameters,\@filtering_param_values,\@filtering_conditions);
combs(0,"",\@NABfiltering_parameters1,\@NABfiltering_param_values1,\@NABfiltering_conditions1);
combs(0,"",\@NABfiltering_parameters2,\@NABfiltering_param_values2,\@NABfiltering_conditions2);

#print("DEBUG: @exe_parameters, @exe_param_values");
#print("DEBUG: @exe_conditions,@filtering_conditions,@NABfiltering_conditions1,@NABfiltering_conditions2\n");

if ($n_cores>1)
{
    $parallel = Parallel::Loops->new($n_cores);
    $parallel->share(\%results);
}


## NAB filtering and parsing
##########################################################################################

if ( ! -f "NAB.vcf")
{
    die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the variant calling step has finished succesfully and try to execute this script again\n";
}

my @NAB_hash_pointers;

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

my @NAB_hash_pointers1;
my @NAB_hash_pointers2;

for (my $i=0; $i<scalar @NABfiltering_conditions1; ++$i)
{
    $NAB_hash_pointers1[$i]=parse_vcf("NAB$sep_param$NABfiltering_conditions1[$i].vcf");
}

for (my $i=0; $i<scalar @NABfiltering_conditions2; ++$i)
{
    $NAB_hash_pointers2[$i]=parse_vcf("NAB$sep_param$NABfiltering_conditions2[$i].vcf");
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

## Old germline filtering
#########################################################################################

#my $filtering_command="$vcf_filt_exe ";
#$filtering_command.="--isvar 1";
#if (!-f "NAB${sep_param}germline.vcf")
#{
#	print("Filtering NAB.vcf to generate NAB${sep_param}germline.vcf\n");
#	## Filter the right vcf file generated ini the outside loop	
#	system("$filtering_command -i NAB.vcf -o NAB${sep_param}germline.vcf");
#	#print("DEBUG: $filtering_command -i A$sep_param$exe_condition.vcf -o A$sep_param$condition.vcf \n");
#}
#else
#{
#	print("NAB${sep_param}germline.vcf had been previously generated and it will be recycled\n");
#}


##my @temp=parse_vcf_name("NAB${sep_param}germline.vcf");
##my %NAB=%{$temp[0]};
##my $nameN=$temp[1];

foreach my $exe_condition (@exe_conditions) ##Options that require to call variants again 
{
#print("DEBUG: Exe condition loop $exe_condition\n");
    if ((!-f "A$sep_param$exe_condition.vcf") || (!-f "B$sep_param$exe_condition.vcf")) ##VCF file we will use for filtering, if it does not exist we have to perform the variant calling
    {
        die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the two scripts are using the same filtering options\n";
    }
    my @current_conditions=@filtering_conditions;
    for (my $i=0; $i<scalar @filtering_conditions;++$i)
    {
        $current_conditions[$i]=[$exe_condition,$filtering_conditions[$i]];
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
print($OFILE "Sample,Condition,A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean\n");
#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean
my $sample=$output_file;
$sample=basename($sample);
$sample=~s/(.*)\..*/$1/;
foreach my $condition (keys %results)
{
    print($OFILE "$sample$OFS$condition$OFS",array_to_string(@nofilt_results),$OFS,array_to_string(@{$results{$condition}}),"\n");
}
close($OFILE);

## MAIN BODY AS FUNCTION FOR PARALLELISM
########################################
sub filter
{
    my $exe_condition;
    my $filtering_condition;
    if (defined $_)
    {
        ($exe_condition,$filtering_condition)=@{$_};
    }
    else
    {
        ($exe_condition,$filtering_condition)=@{$_[0]};
    }
    my $condition="$exe_condition$sep_param$filtering_condition";
    my $filtering_command="$vcf_filt_exe ";
    $filtering_command.=join(" ",split("$sep_value",join(" ",split("$sep_param",$filtering_condition))));

    if (!-f "A$sep_param$condition.vcf")
    {
        print("Filtering A$sep_param$exe_condition.vcf to generate A$sep_param$condition.vcf\n");
## Filter the right vcf file generated in the outside loop	
        system("$filtering_command -i A$sep_param$exe_condition.vcf -o A$sep_param$condition.vcf"); ## I don't think speed would be an issue, thus I'll call the filtering script every time (instead of writting a package and/or functions here)
#print("DEBUG: $filtering_command -i A$sep_param$exe_condition.vcf -o A$sep_param$condition.vcf \n");

    }
    else
    {
        print("A$sep_param$condition.vcf has been previously generated and it will be recycled\n");
    }

    if (!-f "B$sep_param$condition.vcf")
    {
        print("Filtering B$sep_param$exe_condition.vcf to generate B$sep_param$condition.vcf\n");
## Filter the right vcf file generated in the outside loop		
        system("$filtering_command -i B$sep_param$exe_condition.vcf -o B$sep_param$condition.vcf"); ## I don't think speed would be an issue, thus I'll call the filtering script every time (instead of writting a package and/or functions here)
#print("DEBUG: $filtering_command -i B$sep_param$exe_condition.vcf -o B$sep_param$condition.vcf \n");
    }
    else
    {
        print("B$sep_param$condition.vcf has been previously generated and it will be recycled\n");
    }

    my %Afilt=%{parse_vcf("A$sep_param$condition.vcf")};
    my %Bfilt=%{parse_vcf("B$sep_param$condition.vcf")};

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

#Output of list of variants and/or intermediate vcf files
    if($output_list)
    {
        if($output_comprehensive)
        {
            write_variant_list($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.list");
            write_variant_list($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.list");
            write_variant_list($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.list");
            write_variant_list($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.list");
            write_variant_list($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.list");
            write_variant_list($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.list");

        }

        write_variant_list($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.list");
        write_variant_list($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.list");
        write_variant_list($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.list");
        write_variant_list($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.list");

    }

    if($output_vcfs)
    {
        if($output_comprehensive)
        {
            write_variant_vcf($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.vcf","A.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.vcf","B.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.vcf","A.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.vcf","B.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.vcf","A.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.vcf","B.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.vcf","A.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_vcf($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.vcf","B.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.vcf","A.vcf","B.vcf","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.vcf","A.vcf","B.vcf","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.vcf","A.vcf","B.vcf","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
            write_variant_2vcf($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.vcf","A.vcf","B.vcf","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");

        }
        write_variant_2vcf($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.vcf","A.vcf","B.vcf","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.vcf","A.vcf","B.vcf","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.vcf","A.vcf","B.vcf","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");
        write_variant_2vcf($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.vcf","A.vcf","B.vcf","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}");

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

#Mean stats filter NAB		
        my @statsfiltNABU;
        my @statsfiltNABI;

        my ($ref_common_variantsfiltNABU,$ref_different_variantsfiltNABU)=vcf_prune($ref_common_variantsfiltNU,$ref_different_variantsfiltNU,$NAB,\@statsfiltNABU);
        my ($ref_common_variantsfiltNABI,$ref_different_variantsfiltNABI)=vcf_prune($ref_common_variantsfiltNI,$ref_different_variantsfiltNI,$NAB,\@statsfiltNABI);

        my @statsfiltNABmean=(($statsAfiltNAB[0]+$statsBfiltNAB[0])/2.0,($statsAfiltNAB[1]+$statsBfiltNAB[1])/2.0,($statsAfiltNAB[2]+$statsBfiltNAB[2])/2.0);

#Output of list of variants and/or intermediate vcf files
        if($output_list)
        {
            if($output_comprehensive)
            {
                write_variant_list($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
                write_variant_list($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
                write_variant_list($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");
                write_variant_list($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list");

            }

            write_variant_list($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
            write_variant_list($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list");
        }
        if($output_vcfs)
        {
            if($output_comprehensive)
            {
                write_variant_vcf($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","A.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","B.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","A.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_vcf($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","B.vcf","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_2vcf($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","A.vcf","B.vcf","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
                write_variant_2vcf($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","A.vcf","B.vcf","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");


            }
            write_variant_2vcf($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","A.vcf","B.vcf","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");
            write_variant_2vcf($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","A.vcf","B.vcf","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition");

        }

#Store and/or print
        my @statistics=(@statsAfilt,@statsBfilt,@statsfiltU,@statsfiltI,@statsfiltmean,@statsAfiltN,@statsBfiltN,@statsfiltNU,@statsfiltNI,@statsfiltNmean,@statsAfiltNAB,@statsBfiltNAB,@statsfiltNABU,@statsfiltNABI,@statsfiltNABmean);
#Condition,A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean
        $results{"$condition${sep_param}NAB$sep_param$NAB_condition"}=\@statistics;
#print("DEBUG:$condition$OFS",array_to_string(@statistics),"\n");
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
    open(my $OFILE, ">$filename");
    open(my $IFILE, "$vcf");
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

	if (!-f "NAB$sep_param$condition.vcf")
	{
		#print("Filtering NAB$sep_param$exe_condition.vcf to generate NAB$sep_param$condition.vcf\n");
		print("Filtering NAB.vcf to generate NAB$sep_param$condition.vcf\n");
		## Filter the right vcf file generated in the outside loop	
		#system("$filtering_command -i NAB$sep_param$exe_condition.vcf -o NAB$sep_param$condition.vcf");
		system("$filtering_command -i NAB.vcf -o NAB$sep_param$condition.vcf");
		#print("DEBUG: $filtering_command -i A$sep_param$exe_condition.vcf -o A$sep_param$condition.vcf \n");
	}
	else
	{
		print("NAB$sep_param$condition.vcf has been previously generated and it will be recycled\n");
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
            $temp[$j]=~s/0_/${pos_N}_/;
        }
        $array[$i]=\@temp;
    #    print(" after modification ",join(",",@{$array[$i]}),"\n");
    }
    return @array;
   
}
