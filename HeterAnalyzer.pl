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

my $sep_param="#";
my $sep_value=":";
my $OFS=",";
my $FS=",";
my $output_vcf=1;
my $output_list=1;
my $output_comprehensive=1;
my $n_cores=1;
my $covB_sh="covB_half.sh";
my $covN_sh="covN.sh";
#my $covB_sh="covBdoublecheck.sh";

######################################################

##Global Variables
######################################################
my %dictrealnamelist;
my @listnames;
my %dictrealnamevcf;
my @vcfnames;
my @namedvcffiles;

##IO Variables
######################################################
my $execond_inputfile="";
my $filtercond_inputfile="";
my $NABfiltercond_inputfile1="";
my $NABfiltercond_inputfile2="";
my $covBfiltercond_inputfile="";
my $output_file="";
my $SCRIPTSVCF_DIR=$ENV{'SCRIPTSVCF_DIR'};

#Flags
my $help;
my $usage="Usage: $0 [options] -o output_file \n\nWARNING: This is a secondary script that is not inteded to be executed directly.\n\n\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the NAB sample\n\t--NABfilt_cond_inputfile2 : secondary input file for a secondary filter (OR filter implemented in a dirty way) of the NAB sample\n\t--covaltB_cond_inputfile : input file for the filtering taking into account characteristics of the unfiltered in the comparison\n\t--output_dir : output directory for vcf files\n\t--n_cores : number of cores to execute some steps in parallel (requires the perl package Parallel::Loops)\n\t--output_vcf: (bool) generate resulting vcf files or not\n--output_list: (bool) generate resulting list of variants or not\n--comp: (int) indicating the comprehensiveness of the output, 0=no files, 1=only needed files to call variants, 2= all intermediate variants\n\n";
######################################################

######################################################
##MAIN
######################################################
#TODO: check if there is a dictionary with files, read it and rename all files to number.vcf_bkp. Then, when I am going to filter I check if there is the need if the file is there with the _bkp, otherwise overwrite. If there is no need to recalculate, mv to the new name. I may want to rewrite the dict file every time I write a vcf file, to make sure that if I re-start a lost job everything is working fine
##Getopt
######################################################
(! GetOptions(
    'exec_cond_inputfile|e=s' => \$execond_inputfile,
	'filt_cond_inputfile|f=s' => \$filtercond_inputfile,
	'NABfilt_cond_inputfile=s' => \$NABfiltercond_inputfile1,
    'NABfilt_cond_inputfile2=s' => \$NABfiltercond_inputfile2,
    'covaltB_cond_inputfile=s' => \$covBfiltercond_inputfile,
	'output_file|o=s' => \$output_file,
    'n_cores=i' => \$n_cores,
    'output_vcf=i'=>\$output_vcf,
    'output_list=i' => \$output_list,
    'comp=i' => \$output_comprehensive,
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
my @covBfiltering_parameters=("");

my @exe_param_values=([("")]);
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

push(@namedvcffiles,("A.vcf","B.vcf","N.vcf"));

if ( -f "N.vcf")
{
    ($ref_n,$nameN)=parse_vcf_name("N.vcf");
    #print("DEBUG: N genotype name $nameN\n");
    %N=%{$ref_n};
}
else
{
    die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the variant calling step has finished succesfully and try to execute this script again\n";
}

## Main conditions loop
#######################

my $name_condition;
my @filtering_conditions;
my @NABfiltering_conditions1; ##In order to be easily accesed from the filtering functions
my @NABfiltering_conditions2;
my @NABfiltering_conditions;
my @covBfiltering_conditions;
my @exe_conditions;
my %results;
my $bamfiles;
my $parallel;
my %constExeCondContentRefs;
my %nofiltResultsRef;

##Generate all combinations of proposed values for execution and filtering parameters
combs(0,"",\@exe_parameters,\@exe_param_values,\@exe_conditions);
combs(0,"",\@filtering_parameters,\@filtering_param_values,\@filtering_conditions);
combs(0,"",\@NABfiltering_parameters1,\@NABfiltering_param_values1,\@NABfiltering_conditions1);
combs(0,"",\@NABfiltering_parameters2,\@NABfiltering_param_values2,\@NABfiltering_conditions2);
combs(0,"",\@covBfiltering_parameters,\@covBfiltering_param_values,\@covBfiltering_conditions);

#print("DEBUG: @exe_parameters, @exe_param_values");
#print("DEBUG: @exe_conditions,@filtering_conditions,@NABfiltering_conditions1,@NABfiltering_conditions2\n");

##Generating final NAB filtering conditions
###########################################
my $pos=0;
for (my $i=0; $i<scalar @NABfiltering_conditions1; ++$i)
{
    for (my $j=0; $j< scalar @NABfiltering_conditions2; ++$j)
    {
        $NABfiltering_conditions[$pos]="$NABfiltering_conditions1[$i]$sep_param$NABfiltering_conditions2[$j]";
        $pos+=1;
    }
}

if ($n_cores>1)
{
    $parallel = Parallel::Loops->new($n_cores);
    $parallel->share(\%results, \%dictrealnamelist, @listnames, \%dictrealnamevcf, @vcfnames, \@namedvcffiles, \%constExeCondContentRefs, \%nofiltResultsRef);
}

#Parsing covX files
#TODO:I am not 100% sure this is worth it since the forking process takes quite a lot of time to duble-check if I have time
if($n_cores>1 && scalar @exe_conditions > 1)
{
    $parallel->foreach(\@exe_conditions,\&parse_const_execond);
}
else
{
   foreach my $exe_condition (@exe_conditions)
   {
       parse_const_execond($exe_condition);
   }
}

my ($Aexecondname,$Bexecondname);

foreach my $exe_condition (@exe_conditions) ##Options that require to call variants again 
{
    $Aexecondname="A$sep_param$exe_condition.vcf"; ##So far I am not using numbers for execondnames. If I change this, I will need to get the number names here and in the filter function
    $Bexecondname="B$sep_param$exe_condition.vcf";

#print("DEBUG: Exe condition loop $exe_condition\n");
    if ((!-f $Aexecondname) || (!-f $Bexecondname)) ##VCF file we will use for filtering, if it does not exist we have to perform the variant calling
    {
        die "Missing vcf files. Something weird has happend between the execution of the previous script and this one. Check that the two scripts are using the same filtering options\n";
    }
    else
    {
        push(@namedvcffiles,$Aexecondname); ##Fixed name instead of numbers
        push(@namedvcffiles,$Bexecondname);
    }

    my @current_conditions;
    for (my $i=0; $i<scalar @filtering_conditions;++$i)
    {
        $current_conditions[$i]=[$exe_condition,$filtering_conditions[$i]];
    }
#    ##DEBUG
#    foreach my $cond (@current_conditions)
#    {
#        print("DEBUG: ", join(" ",@{$cond}, "\n"));
#    }

    if($n_cores>1 && scalar @current_conditions > 1)
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
$sample=~s/([^.]*)\..*/$1/;
foreach my $condition (keys %results)
{
    print($OFILE "$sample$OFS$condition$OFS",array_to_string(@{$results{$condition}}),"\n");
}
close($OFILE);
my $tag=$exe_conditions[0];
$tag=~s/[^0-9]+//g;
write_vcfrefname_dict("vcfdict.$tag.csv");
write_listrefname_dict("listdict.$tag.csv");
print("Done!\n");

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
   
    my $Aexecondname="A$sep_param$exe_condition.vcf"; ##So far I am not using numbers for execondnames. If I change this I need to obtain the names here from the array (I may make a hash then)
    my $Bexecondname="B$sep_param$exe_condition.vcf";
    my $AcovBname="A$sep_param$exe_condition${sep_param}covBfiltering.tsv"; #Same here, although these are tsvs
    my $BcovBname="B$sep_param$exe_condition${sep_param}covBfiltering.tsv"; 
    my $covNname="covN$sep_param$exe_condition.tsv";

    ##Reading pre-calculated shared data    
    my($refA,$refB)=($constExeCondContentRefs{$Aexecondname}, $constExeCondContentRefs{$Bexecondname});

    my $thisAname=vcf_refname("A$sep_param$condition.vcf",$condition);
    my $thisBname=vcf_refname("B$sep_param$condition.vcf",$condition);
    
    if (!-f $thisAname)
    {
        print("Filtering $Aexecondname to generate $thisAname\n");
        ## Filter the right vcf file generated in the outside loop	
        system("$filtering_command -i $Aexecondname -o $thisAname"); ## I don't think speed would be an issue, thus I'll call the filtering script every time (instead of writting a package and/or functions here)
        #print("DEBUG: $filtering_command -i $Aexecondname -o $thisAname\n");

    }
    else
    {
        print("$thisAname has been previously generated and it will be recycled. WARNING: this may generate inconsistencies if the number or order of parameters are changed between executions\n");
    }

#    print("DEBUG: $Aexecondname generated $thisAname using the filter $condition\n");
 
    if (!-f $thisBname)
    {
        print("Filtering $Bexecondname to generate $thisBname\n");
        ## Filter the right vcf file generated in the outside loop		
        system("$filtering_command -i $Bexecondname -o $thisBname"); ## I don't think speed would be an issue, thus I'll call the filtering script every time (instead of writting a package and/or functions here)
        #print("DEBUG: $filtering_command -i $Bexecondname -o $thisBname \n");
    }
    else
    {
        print("$thisBname has been previously generated and it will be recycled. WARNING: this may generate inconsistencies if the number or order of parameters are changed between executions\n");
    }
    
#    print("DEBUG: $Bexecondname generated $thisBname using the filter $condition\n");

    my %Afilt=%{parse_vcf($thisAname)};
    my %Bfilt=%{parse_vcf($thisBname)};
    my $refdataAfiltcovB=$constExeCondContentRefs{$AcovBname};
    my $refdataBfiltcovB=$constExeCondContentRefs{$BcovBname};
    my $refdatacovN=$constExeCondContentRefs{$covNname};

    #print("DEBUG: reference B $refB, reference A $refA, $filtering_condition\n");
    
    #Compare Afilt with B without filter --> Common variants + %
    my @statsAfilt;
    my ($ref_common_variantsAfilt,$ref_different_variantsAfilt)=vcf_compare_parsed($refB,\%Afilt,\@statsAfilt); ##I have to generate two hashes. One with common variants, the other with non common. Thus, the consecutive filter I can do it towards these new (smallest) hashes.
    
    #Compare Bfiltered with A without filter --> Common variants + %
    my @statsBfilt;
    my ($ref_common_variantsBfilt,$ref_different_variantsBfilt)=vcf_compare_parsed($refA,\%Bfilt,\@statsBfilt);
    #print("DEBUG: statsA of $filtering_condition @statsAfilt\n");

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

    #print(join(",","DEBUG nofilt:",@{$nofiltResultsRef{$exe_condition}},"\n"));
    #print(join(",","DEBUG filt and filtN:",@statsAfilt,@statsBfilt,@statsfiltU,@statsfiltI,@statsAfiltN,@statsBfiltN,@statsfiltNU,@statsfiltNI,@statsfiltNmean,"\n")); #DEBUG
    
    #Variables that will be used just in the first iteraction of ncovB and reused for all NABs
    my @statsAfiltNAB;
    my @statsBfiltNAB;
    my @statsfiltNABU;
    my @statsfiltNABI;
    my ($ref_common_variantsAfiltNAB,$ref_different_variantsAfiltNAB);
    my ($ref_common_variantsBfiltNAB,$ref_different_variantsBfiltNAB);
    my ($ref_common_variantsfiltNABU,$ref_different_variantsfiltNABU);
    my ($ref_common_variantsfiltNABI,$ref_different_variantsfiltNABI);
    my @statsfiltNABmean;
    my @statsfiltNABrefs; #NAB is inside the covB loop while they are independent, in order to have just one set of nested loops instead of two single loops and one set of nested loops (since we have, covB, NAB and their interaction, covBNAB). I save here the ref to an array with stats for each NAB condition
 
    for (my $ncovB=0; $ncovB < scalar @covBfiltering_conditions; ++$ncovB)
    {   
        #CovB loop variables
        my $covBfiltering_condition;
    	my $NcovBcondition;        	
        my @statsAfiltNcovB;
        my @statsBfiltNcovB;
        my @statsfiltNcovBU;
        my @statsfiltNcovBI;
        my ($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB);
        my ($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB);
        my ($ref_common_variantsfiltNcovBU,$ref_different_variantsfiltNcovBU);
        my ($ref_common_variantsfiltNcovBI,$ref_different_variantsfiltNcovBI);
        my @statsfiltNcovBmean;

        $covBfiltering_condition=$covBfiltering_conditions[$ncovB];
    	$NcovBcondition="$exe_condition$sep_param$filtering_condition${sep_param}covB$sep_param$covBfiltering_condition";        	
    
        #AfiltNcovB
        ($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB)=vcf_prune_covB($ref_common_variantsAfiltN,$ref_different_variantsAfiltN,$refdataBfiltcovB,\@statsAfiltNcovB,$covBfiltering_condition);
    
        #BfiltNcovB
        ($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB)=vcf_prune_covB($ref_common_variantsBfiltN,$ref_different_variantsBfiltN,$refdataAfiltcovB,\@statsBfiltNcovB,$covBfiltering_condition);
    
        #Union
        ($ref_common_variantsfiltNcovBU,$ref_different_variantsfiltNcovBU)=vcf_unite_parsed($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,\@statsfiltNcovBU);
        
    	#Intersection
        ($ref_common_variantsfiltNcovBI,$ref_different_variantsfiltNcovBI)=vcf_intersect_parsed($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,\@statsfiltNcovBI);
    
        #Means
        @statsfiltNcovBmean=(($statsAfiltNcovB[0]+$statsBfiltNcovB[0])/2.0,($statsAfiltNcovB[1]+$statsBfiltNcovB[1])/2.0,($statsAfiltNcovB[2]+$statsBfiltNcovB[2])/2.0);
    
        #print(join(",","DEBUG filtNcovB:",@statsAfiltNcovB,@statsBfiltNcovB,@statsfiltNcovBU,@statsfiltNcovBI,@statsfiltNcovBmean,"\n")); #DEBUG
    
        #Output of list of variants and/or intermediate vcf files
        if($output_comprehensive>0)
        {
            if ($output_list)
            {
                write_variant_list($ref_common_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_common.list",$condition);
                write_variant_list($ref_different_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_different.list",$condition);
                write_variant_list($ref_different_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_different.list",$condition);
            }
            if ($output_vcf)
            {
                write_variant_2vcf($ref_common_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_common.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                write_variant_vcf($ref_different_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                write_variant_vcf($ref_different_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
            }
            if($output_comprehensive>1)
            {
                if ($output_list)
                {
                    write_variant_list($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_common.list",$condition);
                    write_variant_list($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.list",$condition);
                    write_variant_list($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.list",$condition);
                }
                if ($output_vcf)
                {
                    write_variant_vcf($ref_common_variantsAfilt,"Afilt$sep_param${condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsBfilt,"Bfilt$sep_param${condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsAfiltN,"AfiltN$sep_param${condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsBfiltN,"BfiltN$sep_param${condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsAfiltNcovB,"AfiltNcovB$sep_param${NcovBcondition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_common_variantsBfiltNcovB,"BfiltNcovB$sep_param${NcovBcondition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsAfilt,"Afilt$sep_param${condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsBfilt,"Bfilt$sep_param${condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsAfiltN,"AfiltN$sep_param${condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_vcf($ref_different_variantsBfiltN,"BfiltN$sep_param${condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltU,"filtU$sep_param${condition}_different.vcf","$thisAname","$thisBname","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNU,"filtNU$sep_param${condition}_different.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNcovBU,"filtNcovBU$sep_param${NcovBcondition}_different.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltI,"filtI$sep_param${condition}_different.vcf","$thisAname","$thisBname","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNI,"filtNI$sep_param${condition}_different.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_different_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_different.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltU,"filtU$sep_param${condition}_common.vcf","$thisAname","$thisBname","## U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltNU,"filtNU$sep_param${condition}_common.vcf","$thisAname","$thisBname","## Somatic U(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltI,"filtI$sep_param${condition}_common.vcf","$thisAname","$thisBname","## I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltNI,"filtNI$sep_param${condition}_common.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                    write_variant_2vcf($ref_common_variantsfiltNcovBI,"filtNcovBI$sep_param${NcovBcondition}_common.vcf","$thisAname","$thisBname","## Somatic I(AB) Filtered with vcfFilterTableV1. Condition ${condition}",$condition);
                }
            }
        }
    
        for (my $nNAB=0; $nNAB< scalar @NABfiltering_conditions; ++$nNAB)
        { 
            #covB*NAB loop variables
            my $NAB_condition;
            my $NAB;
            my @statsAfiltcovBNAB;
            my ($ref_common_variantsAfiltcovBNAB,$ref_different_variantsAfiltcovBNAB);
            my @statsBfiltcovBNAB;
            my ($ref_common_variantsBfiltcovBNAB,$ref_different_variantsBfiltcovBNAB);
            my @statsfiltcovBNABU;
            my @statsfiltcovBNABI;
            my ($ref_common_variantsfiltcovBNABU,$ref_different_variantsfiltcovBNABU);
            my ($ref_common_variantsfiltcovBNABI,$ref_different_variantsfiltcovBNABI);
            my @statsfiltcovBNABmean;

            $NAB_condition=$NABfiltering_conditions[$nNAB];
            $NAB=$constExeCondContentRefs{"${exe_condition}_${NAB_condition}"};

            if ($ncovB==0)
            {    
        #Substract NAB from AfiltN. Compare the results to B without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
                ($ref_common_variantsAfiltNAB,$ref_different_variantsAfiltNAB)=vcf_prune($ref_common_variantsAfiltN,$ref_different_variantsAfiltN,$NAB,\@statsAfiltNAB);
        
        #Substract NAB from BfiltN. Compare the results to A without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
                ($ref_common_variantsBfiltNAB,$ref_different_variantsBfiltNAB)=vcf_prune($ref_common_variantsBfiltN,$ref_different_variantsBfiltN,$NAB,\@statsBfiltNAB);
        
        #Mean stats filter NAB		
                ($ref_common_variantsfiltNABU,$ref_different_variantsfiltNABU)=vcf_prune($ref_common_variantsfiltNU,$ref_different_variantsfiltNU,$NAB,\@statsfiltNABU);
                ($ref_common_variantsfiltNABI,$ref_different_variantsfiltNABI)=vcf_prune($ref_common_variantsfiltNI,$ref_different_variantsfiltNI,$NAB,\@statsfiltNABI);
                @statsfiltNABmean=(($statsAfiltNAB[0]+$statsBfiltNAB[0])/2.0,($statsAfiltNAB[1]+$statsBfiltNAB[1])/2.0,($statsAfiltNAB[2]+$statsBfiltNAB[2])/2.0);
                $statsfiltNABrefs[$nNAB]=[@statsAfiltNAB,@statsBfiltNAB,@statsfiltNABU,@statsfiltNABI,@statsfiltNABmean];

                if($output_comprehensive>1)
                {
                    if($output_list)
                    {
                        write_variant_list($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                        write_variant_list($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                        write_variant_list($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                        write_variant_list($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                        write_variant_list($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    }
                    if($output_vcf)
                    {
                        write_variant_vcf($ref_common_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_vcf($ref_common_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_vcf($ref_different_variantsAfiltNAB,"AfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_vcf($ref_different_variantsBfiltNAB,"BfiltNAB$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_2vcf($ref_different_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_2vcf($ref_different_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition); 
        			    write_variant_2vcf($ref_common_variantsfiltNABU,"filtNABU$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                        write_variant_2vcf($ref_common_variantsfiltNABI,"filtNABI$sep_param${condition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    }
                }
            }

            #print(join(",","DEBUG filtNAB:",@{$statsfiltNABrefs[$nNAB]},"\n")); #DEBUG
            
            #Substract NAB from AfiltNcovB. Compare the results to B without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
            ($ref_common_variantsAfiltcovBNAB,$ref_different_variantsAfiltcovBNAB)=vcf_prune($ref_common_variantsAfiltNcovB,$ref_different_variantsAfiltNcovB,$NAB,\@statsAfiltcovBNAB);
    
            #Substract NAB from BfiltNcovB. Compare the results to A without filter --> Common variants + % ###We want to apply filter to NAB and remove only variants that are Alternative for N
            ($ref_common_variantsBfiltcovBNAB,$ref_different_variantsBfiltcovBNAB)=vcf_prune($ref_common_variantsBfiltNcovB,$ref_different_variantsBfiltNcovB,$NAB,\@statsBfiltcovBNAB);
    
            #Mean stats filter covBNAB    
            ($ref_common_variantsfiltcovBNABU,$ref_different_variantsfiltcovBNABU)=vcf_prune($ref_common_variantsfiltNcovBU,$ref_different_variantsfiltNcovBU,$NAB,\@statsfiltcovBNABU);
            ($ref_common_variantsfiltcovBNABI,$ref_different_variantsfiltcovBNABI)=vcf_prune($ref_common_variantsfiltNcovBI,$ref_different_variantsfiltNcovBI,$NAB,\@statsfiltcovBNABI);
    
            @statsfiltcovBNABmean=(($statsAfiltcovBNAB[0]+$statsBfiltcovBNAB[0])/2.0,($statsAfiltcovBNAB[1]+$statsBfiltcovBNAB[1])/2.0,($statsAfiltcovBNAB[2]+$statsBfiltcovBNAB[2])/2.0);
    
            #print(join(",","DEBUG filtcovBNAB:",@statsAfiltcovBNAB,@statsBfiltcovBNAB,@statsfiltcovBNABU,@statsfiltcovBNABI,@statsfiltcovBNABmean,"\n")); #DEBUG

            #Output of list of variants and/or intermediate vcf files
            
            if($output_comprehensive>0)
            {
                if($output_list)
                {
                    write_variant_list($ref_common_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    write_variant_list($ref_different_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                }
                if($output_vcf)
                {
                    write_variant_2vcf($ref_common_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_vcf($ref_different_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_vcf($ref_different_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                }
            }
            if($output_comprehensive>1)
            {
                if($output_list)
                {
    				write_variant_list($ref_common_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    write_variant_list($ref_common_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                    write_variant_list($ref_different_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                    write_variant_list($ref_different_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.list",$condition);
                    write_variant_list($ref_common_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.list",$condition);
                }
                if($output_vcf)
                {
                    write_variant_vcf($ref_common_variantsAfiltcovBNAB,"AfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_vcf($ref_common_variantsBfiltcovBNAB,"BfiltcovBNAB$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_2vcf($ref_different_variantsfiltcovBNABU,"filtcovBNABU$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB U(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_2vcf($ref_different_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_different.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                    write_variant_2vcf($ref_common_variantsfiltcovBNABI,"filtcovBNABI$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}_common.vcf","$thisAname","$thisBname","## Somatic NAB I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition",$condition);
                }
            }

            #Final variables
            my @statistics;
            #Store and/or print
            @statistics=(@{$nofiltResultsRef{$exe_condition}},@statsAfilt,@statsBfilt,@statsfiltU,@statsfiltI,@statsfiltmean,@statsAfiltN,@statsBfiltN,@statsfiltNU,@statsfiltNI,@statsfiltNmean,@statsAfiltNcovB,@statsBfiltNcovB,@statsfiltNcovBU,@statsfiltNcovBI,@statsfiltNcovBmean,@{$statsfiltNABrefs[$nNAB]},@statsAfiltcovBNAB,@statsBfiltcovBNAB,@statsfiltcovBNABU,@statsfiltcovBNABI,@statsfiltcovBNABmean);
    #Condition,#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean
    
    		$results{"$NcovBcondition${sep_param}NAB$sep_param$NAB_condition"}=\@statistics;
    #print("DEBUG:$condition$OFS",array_to_string(@statistics),"\n");
        }
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

    for ($i=0;$i<scalar @tsv1;$i++)
    {
        chomp($tsv1[$i]);
        my @values=split("\t",$tsv1[$i]);
        $key="$values[0]$OFS$values[1]";
        splice(@values,0,2);
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
    my %common_variants=%{$vcf_1};
    my %private_variants=%{$vcf_2};
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

    #print("DEBUG: covB filtering options $covBfiltering_options\n");    

    my @covBfiltering_params=split($sep_param,$covBfiltering_options);
    foreach my $covBoption (@covBfiltering_params)
    {
        my ($param, $value)= split($sep_value,$covBoption);

        #print("DEBUG: covB option $covBoption, param: $param, value: $value\n");
        
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
    #print("DEBUG:Final filtering options, minimum coverage in B: $min_coverage, maximum number of reads supporting the alternative allele in B: $max_alternative, maximum proportion of alternative alleles in B: $max_propalt\n");
    ##Filtering loop
    ###############
    my $remove_var=0;
    my $remove_private=0;

    foreach my $tsv_variant (keys %{$tsv_todelete})
    {

        $totalreads=0;
        $totalaltreads=0;
        $totalaltreads+=$_ for split(",",${${$tsv_todelete}{$tsv_variant}}[3]); ##Oneline for loop
        $totalreads=$totalaltreads+${${$tsv_todelete}{$tsv_variant}}[2];
        $prop_alts=$totalreads != 0 ? ($totalaltreads+0.0)/$totalreads : 0 ;

        if($totalreads < $min_coverage)
        {
            $remove_var=1;
        }
        else
        {
            $remove_var=0;
            if($totalaltreads>$max_alternative || $prop_alts>$max_propalt)
            {
               $remove_private=1; 
            }
            else
            {
                $remove_private=0;
            }
        }

#DEBUG        
#        my $common=exists($common_variants{$tsv_variant})?"yes":"no";
#        my $private=exists($private_variants{$tsv_variant})?"yes":"no";
#        if(exists($common_variants{$tsv_variant}) || exists($private_variants{$tsv_variant}))
#        {
#            #print("DEBUG: Variant $tsv_variant\n");
#            #print("DEBUG: removevar $remove_var, remove_private $remove_private, common $common, private $private, tag1 $tag1, tag2 $tag2, totalreads $totalreads, totalalternative $totalaltreads, proportion alternative $prop_alts,var 1 ${${$tsv_todelete}{$tsv_variant}}[2], var 2 ${${$tsv_todelete}{$tsv_variant}}[3]\n");###$tsv_todelete[2] Number of bases in the reference allele $tsv_todelete[3]
#        }
#DEBUG
        if ($tag1==0 and $remove_var==1 and exists($common_variants{$tsv_variant}))
        {
            delete($common_variants{$tsv_variant});
            #print("\tDEBUG: deleted in common mutation\n");
        }
        if ($tag2==0 and ($remove_var==1 or $remove_private==1) and exists($private_variants{$tsv_variant}))
        {
            delete($private_variants{$tsv_variant});
            #print("\tDEBUG: deleted in private mutation\n");
        }
        if($tag1==0 and scalar(keys %common_variants)==0)
        {
            $tag1=1;
            #print("\tDEBUG: No more common variants\n");
        }
        if($tag2==0 and scalar(keys %private_variants)==0)
        {
            $tag2=1;
            #print("\tDEBUG: No more private variants\n");
        }
		
		if($tag1==1 and $tag2==1)
        {
            #print("\tDEBUG: No more variants\n");
            last;
        }

    }
    my $n_selvariants=scalar keys %common_variants;
    my $n_filtvariants=scalar keys %private_variants;
    if($n_selvariants+$n_filtvariants==0)
    {
         @{ $ref_statistics }=(0,$n_selvariants,$n_selvariants+$n_filtvariants);
    }
    else
    {
        @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    }
    return (\%common_variants,\%private_variants);
}

#Filter out variants in one hash depending on certain filters from two others and generate statistics
#####################################################################################################
sub vcf_covN_filter
{
    
    my ($ref_datacovN,$covNfiltering_options)=@_;
    my %variants=%{$ref_datacovN};
    my $totalreads=0;
    my $totalaltreads=0;
    my $prop_alts=0;
    
    #Parser for covNfiltering_options
    #################################################################################

    my $min_coverage=0;
    my $max_alternative=9**9**9; ##pseudo +inf
    my $max_propalt=1;

    #print("DEBUG: covN filtering options $covNfiltering_options\n");    

    my @covNfiltering_params=split($sep_param,$covNfiltering_options);
    foreach my $covNoption (@covNfiltering_params)
    {
        my ($param, $value)= split($sep_value,$covNoption);
        $value=~s/0_//; ##They start with 0_ for backwards compatiblity

        #print("DEBUG: covN option $covNoption, param: $param, value: $value\n");
        
        if ($param =~ /--min_coverage/i || $param =~ /--max_coverage/i) ##The second term is included for backwards compatibility and it is misleading
        {
           $min_coverage=$value;
        }
        elsif ($param =~ /--max_alternative/i || $param =~ /--min_reads_alternate/i)
        {
           $max_alternative=$value;
        }
        elsif ($param =~ /--max_propalt/i || $param =~ /--min_freq_alt/i)
        {
           $max_propalt=$value;
        }
        else
        {
            die "The param $param has not been recognized properly in the covN filtering step\n";
        }
    }
    #print("DEBUG:Final filtering options, minimum coverage in B: $min_coverage, maximum number of reads supporting the alternative allele in B: $max_alternative, maximum proportion of alternative alleles in B: $max_propalt, variants with data in covN" . scalar (keys %variants) . "\n");
    
    ##Filtering loop
    ###############

    foreach my $variant (keys %variants)
    {

        $totalreads=0;
        $totalaltreads=0;
        $totalaltreads+=$_ for split(",",${$variants{$variant}}[3]); ##Oneline for loop
        $totalreads=$totalaltreads+${$variants{$variant}}[2];
        $prop_alts=$totalreads != 0 ? ($totalaltreads+0.0)/$totalreads : 0 ;

       # print("DEBUG: Variant $variant\n");
       # print("DEBUG: totalreads $totalreads, totalalternative $totalaltreads, proportion alternative $prop_alts\n");

        if($totalreads < $min_coverage || ($totalaltreads>$max_alternative && $prop_alts>$max_propalt))
        {
           # print("\tDEBUG: covN variant will be kept, to be removed from A or B if present\n");
        }
        else
        { 
            #print("\tDEBUG: good N variant, removed from the list\n");
            delete($variants{$variant}); 
        }
    }
    return(\%variants);
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
    my ($ref_hash,$filename,$filter)=@_;
    my $FILE=openwrite_list_refname("$filename",$filter);
    (! defined $FILE) and return;
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
    my ($ref_hash,$filename,$vcf,$comment,$filter)=@_;
    my $OFILE=openwrite_vcf_refname("$filename",$filter);
    (! defined $OFILE) and return;
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
    my ($ref_hash,$filename,$vcf,$vcf2,$comment,$filter)=@_;
    my $OFILE=openwrite_vcf_refname("$filename",$filter);
    (! defined $OFILE) and return;
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

    my $thisNABname="NAB$sep_param$condition.vcf"; #If I ever include this in the dict, I have to remove the next line
    push(@namedvcffiles,$thisNABname);
	
    if (!-f $thisNABname)
	{
		#print("Filtering NAB$sep_param$exe_condition.vcf to generate $thisNABname\n");
		print("Filtering NAB.vcf to generate $thisNABname\n");
		## Filter the right vcf file generated in the outside loop	
		#system("$filtering_command -i NAB$sep_param$exe_condition.vcf -o NAB$sep_param$condition.vcf");
		system("$filtering_command -i NAB.vcf -o $thisNABname");
	}
	else
	{
		print("$thisNABname has been previously generated and it will be recycled\n");
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

# Opens a vcf file to write, using a number instead of its original name
# It calculates a tag eliminating all non-numeric characters from the filter name
# This is necessary for the parallel version since processes cannot share variables (they are only shared at the end)
# A mapping between numbers and original names is printed by write_vcfrefname_dict
# If the file is already in the hash/array system, returns undef
# ######################################################################

sub openwrite_vcf_refname
{
    my ($filename,$filter)=@_;
    $filter=~s/[^0-9]//g;
    if(exists $dictrealnamevcf{$filename} && $dictrealnamevcf{$filename}=~/^$filter/)
    {
        warn "DEBUG: skipping file $filename since it had been previously written by this process\n";
        return undef;
    }
    my $n_name=push (@vcfnames,$filename);
    my $realfilename=$filter . $n_name . ".vcf";
    $dictrealnamevcf{$filename}=$realfilename;
    open (my $FILE, ">$realfilename") or die "Error opening $realfilename, in place of $filename\n";
    return $FILE;
}

# Returns the new vcf filename to write, using a number instead of its original name
# It calculates a tag eliminating all non-numeric characters from the filter name
# This is necessary for the parallel version since processes cannot share variables (they are only shared at the end)
# A mapping between numbers and original names has to be then printed using write_vcfrefname_dict
# If the file is already in the hash/array system, returns the name previously stored
# ###############################################################################################

sub vcf_refname
{
    my ($filename,$filter)=@_;
    $filter=~s/[^0-9]//g;
    if(exists $dictrealnamevcf{$filename} && $dictrealnamevcf{$filename}=~/^$filter/)
    {
        warn "DEBUG: this file, $filename had been previously written by this process\n";
        return $dictrealnamevcf{$filename}
    }
    my $n_name=push (@vcfnames,$filename);
    my $realfilename=$filter . $n_name . ".vcf";
    $dictrealnamevcf{$filename}=$realfilename; 
    return $realfilename;
}

# Opens a list file to write, using a number instead of its original name
# It calculates a tag eliminating all non-numeric characters from the filter name
# This is necessary for the parallel version since processes cannot share variables (they are only shared at the end)
# A mapping between numbers and original names is printed by write_listrefname_dict
# ######################################################################
sub openwrite_list_refname
{
    my ($filename,$filter)=@_;
    $filter=~s/[^0-9]//g;
    if(exists $dictrealnamelist{$filename} && $dictrealnamelist{$filename}=~/^$filter/)
    {
        warn "DEBUG: skipping file $filename since it had been previously written by this process\n";
        return undef; 
    }
    my $n_name=push(@listnames,$filename);
    my $realfilename=$filter . $n_name . ".list";
    $dictrealnamelist{$filename}=$realfilename;
    open (my $FILE, ">$realfilename") or die "Error opening $realfilename, in place of $filename\n";
    return $FILE;
}


# Generates the output file with the mapping of vcf file numbers and ideal names
# ##############################################################################
sub write_vcfrefname_dict
{
    my ($filename)=@_;
    open (my $FILE, ">$filename");
    foreach my $filename (keys %dictrealnamevcf)
    {
        print($FILE "$filename$OFS$dictrealnamevcf{$filename}\n");
    }
    foreach my $file (@namedvcffiles)
    {
        print($FILE "$file$OFS$file\n");
    }
    close($FILE);
}

# Generates the output file with the mapping of list file numbers and ideal names
# ##############################################################################

sub write_listrefname_dict
{
    my ($filename)=@_;
    open (my $FILE, ">$filename");
    foreach my $filename (keys %dictrealnamelist)
    {
        print($FILE "$filename$OFS$dictrealnamelist{$filename}\n");
    }
    close($FILE);
}

# Parses all const_execond
# ###################################################################

sub parse_const_execond
{
    my $exe_condition;
    if (defined $_)
    {
        $exe_condition=$_;
    }
    else
    {
        $exe_condition=$_[0];
    }
    
    my $Aexecondname="A$sep_param$exe_condition.vcf";
    my $Bexecondname="B$sep_param$exe_condition.vcf";
    my $AcovBname="A$sep_param$exe_condition${sep_param}covBfiltering.tsv";
    my $BcovBname="B$sep_param$exe_condition${sep_param}covBfiltering.tsv";
    my $covNname="covN$sep_param$exe_condition.tsv";
    
    (! -f $covNname) and die "Error calculating coverage of N";
    (! -f $AcovBname ) or (! -f $BcovBname) and die "Error calculating coverages of A and B";
    (! -f $Aexecondname) or  (! -f $Bexecondname) and die "Missing vcf files for constant files. Check the variant calling step";

    $constExeCondContentRefs{$AcovBname}=parse_tsv($AcovBname);
    $constExeCondContentRefs{$BcovBname}=parse_tsv($BcovBname);
    $constExeCondContentRefs{$covNname}=parse_tsv($covNname);
    
    #TODO: I am not parallelizing this, but I am getting a lot of benefits from reusing it. I think that parallalelizing it would spend more extra time spliting the fork than the actual benefit due to the several processes
    foreach my $cond (@NABfiltering_conditions)
    {
        $constExeCondContentRefs{"${exe_condition}_${cond}"}=vcf_covN_filter($constExeCondContentRefs{$covNname},$cond);
    }
    
    $constExeCondContentRefs{$Aexecondname}=parse_vcf($Aexecondname);
    $constExeCondContentRefs{$Bexecondname}=parse_vcf($Bexecondname);
#    $constExeCondContentRefs{"${Aexecondname}_N"}=vcf_prune_single($constExeCondContentRefs{$Aexecondname},\%N); #Not needed
#    $constExeCondContentRefs{"${Bexecondname}_N"}=vcf_prune_single($constExeCondContentRefs{$Bexecondname},\%N);
    $nofiltResultsRef{$exe_condition}=[scalar keys %{$constExeCondContentRefs{$AcovBname}}, scalar keys %{$constExeCondContentRefs{$BcovBname}}, scalar keys %N, scalar keys %{vcf_prune_single($constExeCondContentRefs{$Aexecondname},\%N)}, scalar keys %{vcf_prune_single($constExeCondContentRefs{$Bexecondname},\%N)}];
}
