#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;
use Env;
use Sort::Key::Natural;
use Sort::Key::Maker nat_i_sorter => qw(nat integer);
use Bio::DB::HTS::Tabix;
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
my $popAFfiltercond_inputfile="";
my $output_file="";
my $SCRIPTSVCF_DIR=$ENV{'SCRIPTSVCF_DIR'};
my $GNOMAD=$ENV{'GNOMAD'};

#Flags
my $help;
my $usage="Usage: $0 [options] -o output_file \n\nWARNING: This is a secondary script that is not inteded to be executed directly.\n\n\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the NAB sample\n\t--NABfilt_cond_inputfile2 : secondary input file for a secondary filter (OR filter implemented in a dirty way) of the NAB sample\n\t--covaltB_cond_inputfile : input file for the filtering taking into account characteristics of the unfiltered in the comparison\n\t--popAF_cond_inputfile: input file for the filter of population allele frequencies using gnomAD\n\t--output_dir : output directory for vcf files\n\t--n_cores : number of cores to execute some steps in parallel (requires the perl package Parallel::Loops)\n\t--output_vcf: (bool) generate resulting vcf files or not\n--output_list: (bool) generate resulting list of variants or not\n--comp: (int) indicating the comprehensiveness of the output, 0=no files, 1=only needed files to call variants, 2= all intermediate variants\n\n";
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
    'popAF_cond_inputfile=s' => \$popAFfiltercond_inputfile,
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
my @popAFfiltering_parameters=("");

my @exe_param_values=([("")]);
my @filtering_param_values=([("")]);
my @NABfiltering_param_values1=([("")]);
my @NABfiltering_param_values2=([("")]);
my @covBfiltering_param_values=([("")]);
my @popAFfiltering_param_values=([("")]);

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

if ($popAFfiltercond_inputfile ne "")
{
    @popAFfiltering_parameters=();
    parse_parameters_values($popAFfiltercond_inputfile,\@popAFfiltering_parameters,\@popAFfiltering_param_values);
}

my $vcf_filt_exe="$SCRIPTSVCF_DIR/vcf_filtering.pl";

if (! -f "$vcf_filt_exe")
{
    die "The executable vcf_filtering.pl can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR indicates the directory with this package of scripts\n";
}

unless (-s "$GNOMAD" && -s "${GNOMAD}.tbi")
{
    die "The gnomAD data can't be found at $GNOMAD and $GNOMAD.tbi\n";
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
my @popAFfiltering_conditions;
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
combs(0,"",\@popAFfiltering_parameters,\@popAFfiltering_param_values,\@popAFfiltering_conditions);

##Generate headers of statistics that change in number and name depending on the input filtering conditions
my @HeaderStatsfiltNPAF,
my @HeaderStatsfiltNcovBPAF;
my @HeaderStatsfiltNABPAF;
my @HeaderStatsfiltNABcovBPAF;
my @filtNPAFnames=qw(AfiltN_Priv_PAF_ BfiltN_Priv_PAF_ filtN_U_PAF_);
my @filtNcovBPAFnames=qw(AfiltNcovB_Priv_PAF_ BfiltNcovB_Priv_PAF_ filtNcovB_U_PAF_);
my @filtNABPAFnames=qw(AfiltNAB_Priv_PAF_ BfiltNAB_Priv_PAF_ filtNAB_U_PAF_);
my @filtNABcovBPAFnames=qw(AfiltNABcovB_Priv_PAF_ BfiltNABcovB_Priv_PAF_ filtNABcovB_U_PAF_);

@HeaderStatsfiltNPAF=makeHeaderPAF(\@filtNPAFnames,\@popAFfiltering_conditions);
@HeaderStatsfiltNcovBPAF=makeHeaderPAF(\@filtNcovBPAFnames,\@popAFfiltering_conditions);
@HeaderStatsfiltNABPAF=makeHeaderPAF(\@filtNABPAFnames,\@popAFfiltering_conditions);
@HeaderStatsfiltNABcovBPAF=makeHeaderPAF(\@filtNABcovBPAFnames,\@popAFfiltering_conditions);

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

#Parsing covX files and population allele frequency data for variants of each exe_condition
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

##Making a block here to use the no warnings qw. It complains because of the # but they are valid
{
    no warnings 'qw';
    print($OFILE join(",",qw(Sample Condition A_# B_# N_# AN_# BN_# Afilt_prop Afilt_N Afilt_# Bfilt_prop Bfilt_N Bfilt_# filt_propU filt_NU filt_#U filt_propI filt_NI filt_#I filt_prop_mean filt_N_mean filt_#_mean AfiltN_prop AfiltN_N AfiltN_# BfiltN_prop BfiltN_N BfiltN_# filtN_propU filtN_NU filtN_#U filtN_propI filtN_NI filtN_#I filtN_prop_mean filtN_N_mean filtN_#_mean),@HeaderStatsfiltNPAF,qw(AfiltNcovB_prop AfiltNcovB_N AfiltNcovB_# BfiltNcovB_prop BfiltNcovB_N BfiltNcovB_# filtNcovB_propU filtNcovB_NU filtNcovB_#U filtNcovB_propI filtNcovB_NI filtNcovB_#I filtNcovB_prop_mean filtNcovB_N_mean filtNcovB_#_mean),@HeaderStatsfiltNcovBPAF,qw(AfiltNAB_prop AfiltNAB_N AfiltNAB_# BfiltNAB_prop BfiltNAB_N BfiltNAB_# filtNAB_propU filtNAB_NU filtNAB_#U filtNAB_propI filtNAB_NI filtNAB_#I filtNAB_prop_mean filtNAB_N_mean filtNAB_#_mean),@HeaderStatsfiltNABPAF,qw(AfiltNABcovB_prop AfiltNABcovB_N AfiltNABcovB_# BfiltNABcovB_prop BfiltNABcovB_N BfiltNABcovB_# filtNABcovB_propU filtNABcovB_NU filtNABcovB_#U filtNABcovB_propI filtNABcovB_NI filtNABcovB_#I filtNABcovB_prop_mean filtNABcovB_N_mean filtNABcovB_#_mean),@HeaderStatsfiltNABcovBPAF,qw(AfiltNABcovBPAF_prop AfiltNABcovBPAF_N AfiltNABcovBPAF_# BfiltNABcovBPAF_prop BfiltNABcovBPAF_N BfiltNABcovBPAF_# filtNABcovBPAF_propU filtNABcovBPAF_NU filtNABcovBPAF_#U filtNABcovBPAF_propI filtNABcovBPAF_NI filtNABcovBPAF_#I filtNABcovBPAF_prop_mean filtNABcovBPAF_N_mean filtNABcovBPAF_#_mean),"\n"));
}
    #Condition,#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,(nPAFfreqs) x (AfiltN_Priv_PAF_XX,BfiltN_Priv_PAF_XX,filtN_U_PAF_XX),AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,(nPAFfreqs) x (AfiltNcovB_Priv_PAF_XX,BfiltNcovB_Priv_PAF_XX,filtNcovB_U_PAF_XX),AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,(nPAFfreqs) x (AfiltNAB_Priv_PAF_XX,BfiltNAB_Priv_PAF_XX,filtNAB_U_PAF_XX),AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean,(nPAFfreqs) x (AfiltNABcovB_Priv_PAF_XX,BfiltNABcovB_Priv_PAF_XX,filtNABcovB_U_PAF_XX),AfiltNABcovBPAF_prop,AfiltNABcovBPAF_N,AfiltNABcovBPAF_#,BfiltNABcovBPAF_prop,BfiltNABcovBPAF_N,BfiltNABcovBPAF_#,filtNABcovBPAF_propU,filtNABcovBPAF_NU,filtNABcovBPAF_#U,filtNABcovBPAF_propI,filtNABcovBPAF_NI,filtNABcovBPAF_#I,filtNABcovBPAF_prop_mean,filtNABcovBPAF_N_mean,filtNABcovBPAF_#_mean

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

    #First group of PAF statistics, for filtN
    my @statsfiltNPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltN, $ref_different_variantsBfiltN, $ref_common_variantsfiltNU);
    
    #Variables that will be used just in the first iteraction of ncovB and reused for all NABs
    my @statsAfiltNAB;
    my @statsBfiltNAB;
    my @statsfiltNABU;
    my @statsfiltNABI;
    my @statsfiltNABPAF;
    my ($ref_common_variantsAfiltNAB,$ref_different_variantsAfiltNAB);
    my ($ref_common_variantsBfiltNAB,$ref_different_variantsBfiltNAB);
    my ($ref_common_variantsfiltNABU,$ref_different_variantsfiltNABU);
    my ($ref_common_variantsfiltNABI,$ref_different_variantsfiltNABI);
    my @statsfiltNABmean;
    my @statsfiltNABrefs; #NAB is inside the covB loop while they are independent, in order to have just one set of nested loops instead of two single loops and one set of nested loops (since we have, covB, NAB and their interaction, covBNAB). I save here the ref to an array with stats for each NAB condition
    my @statsfiltNABPAFrefs; #Like above
 
    for (my $ncovB=0; $ncovB < scalar @covBfiltering_conditions; ++$ncovB)
    {   
        #CovB loop variables
        my $covBfiltering_condition;
    	my $NcovBcondition;        	
        my @statsAfiltNcovB;
        my @statsBfiltNcovB;
        my @statsfiltNcovBU;
        my @statsfiltNcovBI;
        my @statsfiltNcovBPAF;
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
        
        #PAF
        @statsfiltNcovBPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltNcovB, $ref_different_variantsBfiltNcovB, $ref_common_variantsfiltNcovBU);

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
            my @statsfiltcovBNABPAF;
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
        
                @statsfiltNABPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltNAB, $ref_different_variantsBfiltNAB, $ref_common_variantsfiltNABU);
                $statsfiltNABPAFrefs[$nNAB]=\@statsfiltNABPAF;

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
    
            @statsfiltcovBNABPAF=getPAFStats(\@popAFfiltering_conditions, $constExeCondContentRefs{"${exe_condition}_PAF"},$ref_different_variantsAfiltcovBNAB, $ref_different_variantsBfiltcovBNAB, $ref_common_variantsfiltcovBNABU);

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
            
            for (my $nPAF=0; $nPAF< scalar @popAFfiltering_conditions; ++$nPAF)
            {
                my $PAF_condition;
                
                my @statsAfiltcovBNABPAF;                
                my ($ref_common_variantsAfiltcovBNABPAF,$ref_different_variantsAfiltcovBNABPAF);
                
                my @statsBfiltcovBNABPAF;
                my ($ref_common_variantsBfiltcovBNABPAF,$ref_different_variantsBfiltcovBNABPAF);

                my @statsfiltcovBNABUPAF;
                my @statsfiltcovBNABIPAF;
                my ($ref_common_variantsfiltcovBNABUPAF,$ref_different_variantsfiltcovBNABUPAF);
                my ($ref_common_variantsfiltcovBNABIPAF,$ref_different_variantsfiltcovBNABIPAF);
                my @statsfiltcovBNABmeanPAF;
    
                $PAF_condition=$popAFfiltering_conditions[$nPAF];

                #Substract PAF from AfiltcovBNAB. Compare the results to B without filter --> Common variants + %
                ($ref_common_variantsAfiltcovBNABPAF,$ref_different_variantsAfiltcovBNABPAF)=filter_with_PAF($ref_common_variantsAfiltcovBNAB,$ref_different_variantsAfiltcovBNAB,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsAfiltcovBNABPAF);

                #Substract PAF from BfiltcovBNAB. Compare the results to A without filter --> Common variants + %
                ($ref_common_variantsBfiltcovBNABPAF,$ref_different_variantsBfiltcovBNABPAF)=filter_with_PAF($ref_common_variantsBfiltcovBNAB,$ref_different_variantsBfiltcovBNAB,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsBfiltcovBNABPAF);

        
                #Mean stats filter covBNAB    
                ($ref_common_variantsfiltcovBNABUPAF,$ref_different_variantsfiltcovBNABUPAF)=filter_with_PAF($ref_common_variantsfiltcovBNABU,$ref_different_variantsfiltcovBNABU,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsfiltcovBNABUPAF);
                
                ($ref_common_variantsfiltcovBNABIPAF,$ref_different_variantsfiltcovBNABIPAF)=filter_with_PAF($ref_common_variantsfiltcovBNABI,$ref_different_variantsfiltcovBNABI,$constExeCondContentRefs{"${exe_condition}_${PAF_condition}"},\@statsfiltcovBNABIPAF);
        
                @statsfiltcovBNABmeanPAF=(($statsAfiltcovBNABPAF[0]+$statsBfiltcovBNABPAF[0])/2.0,($statsAfiltcovBNABPAF[1]+$statsBfiltcovBNABPAF[1])/2.0,($statsAfiltcovBNABPAF[2]+$statsBfiltcovBNABPAF[2])/2.0);

#                #Output of list of variants and/or intermediate vcf files
                
                if($output_comprehensive>0)
                {
                    if($output_list)
                    {
                        write_variant_list($ref_common_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                        write_variant_list($ref_different_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                    }
                    if($output_vcf)
                    {
                        write_variant_2vcf($ref_common_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisAname","$thisBname","## Somatic covBNABPAF U(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_vcf($ref_different_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_vcf($ref_different_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                    }
                }
                if($output_comprehensive>1)
                {
                    if($output_list)
                    {
        				write_variant_list($ref_common_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                        write_variant_list($ref_common_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                        write_variant_list($ref_different_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                        write_variant_list($ref_different_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.list",$condition);
                        write_variant_list($ref_common_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.list",$condition);
                    }
                    if($output_vcf)
                    {
                        write_variant_vcf($ref_common_variantsAfiltcovBNABPAF,"AfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisAname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_vcf($ref_common_variantsBfiltcovBNABPAF,"BfiltcovBNABPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisBname","##Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_2vcf($ref_different_variantsfiltcovBNABUPAF,"filtcovBNABUPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisAname","$thisBname","## Somatic covBNABPAF U(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_2vcf($ref_different_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_different.vcf","$thisAname","$thisBname","## Somatic covBNABPAF I(AB) Filtered with vcfFilterTableV1. Condition ${NcovBcondition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                        write_variant_2vcf($ref_common_variantsfiltcovBNABIPAF,"filtcovBNABIPAF$sep_param${NcovBcondition}${sep_param}NAB$sep_param${NAB_condition}${sep_param}PAF$sep_param${PAF_condition}_common.vcf","$thisAname","$thisBname","## Somatic covBNABPAF I(AB) Filtered with vcfFilterTableV1. Condition ${condition}${sep_param}NAB$sep_param$NAB_condition${sep_param}PAF$sep_param${PAF_condition}",$condition);
                    }
                }
                #Final variables
                my @statistics;
                #Store and/or print
                @statistics=(@{$nofiltResultsRef{$exe_condition}},@statsAfilt,@statsBfilt,@statsfiltU,@statsfiltI,@statsfiltmean,@statsAfiltN,@statsBfiltN,@statsfiltNU,@statsfiltNI,@statsfiltNmean,@statsfiltNPAF,@statsAfiltNcovB,@statsBfiltNcovB,@statsfiltNcovBU,@statsfiltNcovBI,@statsfiltNcovBmean,@statsfiltNcovBPAF,@{$statsfiltNABrefs[$nNAB]},@{$statsfiltNABPAFrefs[$nNAB]},@statsAfiltcovBNAB,@statsBfiltcovBNAB,@statsfiltcovBNABU,@statsfiltcovBNABI,@statsfiltcovBNABmean,@statsfiltcovBNABPAF,@statsAfiltcovBNABPAF, @statsBfiltcovBNABPAF, @statsfiltcovBNABU, @statsfiltcovBNABI, @statsfiltcovBNABmean);
        #Condition,#A_#,B_#,N_#,AN_#,BN_#,Afilt_prop,Afilt_N,Afilt_#,Bfilt_prop,Bfilt_N,Bfilt_#,filt_propU,filt_NU,filt_#U,filt_propI,filt_NI,filt_#I,filt_prop_mean,filt_N_mean,filt_#_mean,AfiltN_prop,AfiltN_N,AfiltN_#,BfiltN_prop,BfiltN_N,BfiltN_#,filtN_propU,filtN_NU,filtN_#U,filtN_propI,filtN_NI,filtN_#I,filtN_prop_mean,filtN_N_mean,filtN_#_mean,(nPAFfreqs) x (AfiltN_Priv_PAF_XX,BfiltN_Priv_PAF_XX,filtN_U_PAF_XX),AfiltNcovB_prop,AfiltNcovB_N,AfiltNcovB_#,BfiltNcovB_prop,BfiltNcovB_N,BfiltNcovB_#,filtNcovB_propU,filtNcovB_NU,filtNcovB_#U,filtNcovB_propI,filtNcovB_NI,filtNcovB_#I,filtNcovB_prop_mean,filtNcovB_N_mean,filtNcovB_#_mean,(nPAFfreqs) x (AfiltNcovB_Priv_PAF_XX,BfiltNcovB_Priv_PAF_XX,filtNcovB_U_PAF_XX),AfiltNAB_prop,AfiltNAB_N,AfiltNAB_#,BfiltNAB_prop,BfiltNAB_N,BfiltNAB_#,filtNAB_propU,filtNAB_NU,filtNAB_#U,filtNAB_propI,filtNAB_NI,filtNAB_#I,filtNAB_prop_mean,filtNAB_N_mean,filtNAB_#_mean,(nPAFfreqs) x (AfiltNAB_Priv_PAF_XX,BfiltNAB_Priv_PAF_XX,filtNAB_U_PAF_XX),AfiltNABcovB_prop,AfiltNABcovB_N,AfiltNABcovB_#,BfiltNABcovB_prop,BfiltNABcovB_N,BfiltNABcovB_#,filtNABcovB_propU,filtNABcovB_NU,filtNABcovB_#U,filtNABcovB_propI,filtNABcovB_NI,filtNABcovB_#I,filtNABcovB_prop_mean,filtNABcovB_N_mean,filtNABcovB_#_mean,(nPAFfreqs) x (AfiltNABcovB_Priv_PAF_XX,BfiltNABcovB_Priv_PAF_XX,filtNABcovB_U_PAF_XX),AfiltNABcovBPAF_prop,AfiltNABcovBPAF_N,AfiltNABcovBPAF_#,BfiltNABcovBPAF_prop,BfiltNABcovBPAF_N,BfiltNABcovBPAF_#,filtNABcovBPAF_propU,filtNABcovBPAF_NU,filtNABcovBPAF_#U,filtNABcovBPAF_propI,filtNABcovBPAF_NI,filtNABcovBPAF_#I,filtNABcovBPAF_prop_mean,filtNABcovBPAF_N_mean,filtNABcovBPAF_#_mean
        
        		$results{"$NcovBcondition${sep_param}NAB$sep_param$NAB_condition"}=\@statistics;
        #print("DEBUG:$condition$OFS",array_to_string(@statistics),"\n");
                #HERE
            }##TODO working here
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

#Parse vcf: Parses CHROM$OFSPOS as key and REF$OFSALT as value
# Genotype-aware
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
    my $value; 

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            chomp($vcf1[$i]);
            $key=$value=$vcf1[$i];
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/; ####TODO: This may be quicker with a split+join strategy. To check it
            $value=~s/^[^\t]+\t[^\t]+\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2/;
            $hash{$key}=$value;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return \%hash;

}


#Parse complete vcf: Parses CHROM$OFSPOS as key and the whole entry as value
# Genotype-aware
###########################################
sub parse_vcf_complete
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    my $value; 

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            chomp($vcf1[$i]);
            $key=$value=$vcf1[$i];
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
            $hash{$key}=$value;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return \%hash;

}

#Parse tsv
# WARNING: genotype information in TSV-related structures is different than in vcf-related ones. Here I used an array ref, with a different set of information. The format here is REF ALT1,ALTN READS READSALT1,READSALTN.
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
# Genotype-aware
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
    my $value;
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
            chomp($vcf1[$i]);
            $key=$value=$vcf1[$i];
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/; ####TODO: This may be quicker with a split+join strategy. To check it
            $value=~s/^[^\t]+\t[^\t]+\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2/;
            $hash{$key}=$value;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return (\%hash,$name);

}

##TODO: This could be improved in terms on speed following the strategy of the next subroutine. I'm not using this one any more, so I won't update it... at least so far.
##Compare two vcf files and generate statistics
###############################################
#sub vcf_compare2
#{
#    my ($vcf1_file,$vcf2_file,$ref_statistics)=@_;
#
#    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
#    open(my $VCF2,$vcf2_file) or die "The file $vcf2_file is not located in the output folder. Please, check if the variant caller has successfully finished";
#    my @vcf1=<$VCF1>;
#    my @vcf2=<$VCF2>;
#    close($VCF1);
#    close($VCF2);
#    clean_vcfcontent(\@vcf1);
#    clean_vcfcontent(\@vcf2);
#    my %variants1=%{variants_to_hash(\@vcf1)};
#    my %variants2=%{variants_to_hash(\@vcf2)};
#    my %commonvariants;
#    my $n1=scalar(keys %variants1);
#    my $n2=scalar(keys %variants2);
#
#    if ($n1<$n2)
#    {
#        foreach my $variant (keys %variants1)
#        {
#            if (exists $variants2{$variant})
#            {
#                $commonvariants{$variant}=1;
##print("DEBUG: New common variant $variant");
#            }
#        }
#    }
#    else
#    {
#        foreach my $variant (keys %variants2)
#        {
#            if (exists $variants1{$variant})
#            {
#                $commonvariants{$variant}=1;
##print("DEBUG: New common variant $variant");
#            }
#        }
#    }
#    @{ $ref_statistics }=($n1,$n2,scalar (keys %commonvariants));
#    return \%commonvariants;
#
#}

#Compare two variant hashes, return a hash with commonvariants and another one with variants in the problem
# that aren't present in the reference and generate statistics
# Genotype-aware
##############################################
sub vcf_compare_parsed
{
    my ($vcf_reference,$vcf_problem,$ref_statistics)=@_;
    my %variants2=%{$vcf_problem};
    my %commonvariants;
    my $n1=scalar(keys %{$vcf_reference});
    my $n2=scalar(keys %variants2);
    my @refAlts;
    my @probAlts;
    my @diffprobAlts;
    my @cAlts;
    my $found;

    print("\nDEBUG: vcf_compare_parsed\n\n");
    foreach my $variant (keys %{$vcf_reference})
    {
        if (exists $variants2{$variant})
        {
            if ($vcf_reference->{$variant} eq $variants2{$variant})
            {
                $commonvariants{$variant}=$variants2{$variant};
                delete $variants2{$variant}; ###This is the local copy of the original hash that will be returned
                #print("DEBUG: New common variant $variant");
            }
            else ## We need to assess if some alleles are in common or none of them
            {
                @refAlts=split($OFS,$vcf_reference->{$variant});
                @probAlts=split($OFS, $variants2{$variant});
                if(scalar @refAlts == 2 && scalar @probAlts == 2) ##Same variant but different alleles
                {
                    print("DEBUG: Same variant but different alleles, ".$vcf_reference->{$variant}." vs. $variants2{$variant}\n");
                    continue;
                }
                 
                @cAlts=($probAlts[0]); ##Will contain the set of alternatives that are in common
                @diffprobAlts=($probAlts[0]); ##Will contain the set of alternatives that are only present in the problem (we don't care about the ones only present in the reference, i.e., unfiltered sample
               
                ##We compare the alternative alleles 
                ##These loops start in 1, since the 0 element is the reference allele
                for (my $i=1; $i<scalar @probAlts; ++$i)
                {
                    $found=0;
                    for (my $j=1; $j<scalar @refAlts; ++$j)
                    {
                        if($probAlts[$i] eq $refAlts[$j])##This will generate some extra comparisons if pairs have been already found. However, the number of comparisons is very small, I don't think generating extra hashes will save much time
                        {
                            push(@cAlts,$probAlts[$i]);
                            $found=1;
                            last;
                        }
                        
                    }
                    unless($found==1)
                    {
                        push(@diffprobAlts,$probAlts[$i]);
                    }
                }
                ##If there are alleles in common, we add this variant to the list of common variants, but only with the relevant allele
                if(scalar @cAlts >1) ##The element 0 is the reference again
                {
                    $commonvariants{$variant}=join($OFS,@cAlts);
                    print("DEBUG: there were a number of common alleles, $commonvariants{$variant}\n");
                }
                ##If there are alleles that were in the problem set and are not in common, we update them in the problem set
                if(scalar @diffprobAlts >1)
                {
                    $variants2{$variant}=join($OFS,@diffprobAlts);
                    print("DEBUG: there were a number of alleles that were not in the reference, $variants2{$variant}\n");
                }
                else ##There were extra alleles in the reference, but not in the problem. Thus, we can eliminate this variant from the problem
                {
                    print("DEBUG: there were extra alleles in the reference, but not in the problem\n");
                    delete $variants2{$variant};
                }
                
            }

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
# Genotype-aware
########################################################################
sub vcf_prune
{
    my ($vcf_1,$vcf_2,$vcf_todelete,$ref_statistics)=@_;
    my %variants1=%{$vcf_1};
    my %variants2=%{$vcf_2};
    my $tag1=0;
    my $tag2=0;
    my @refAlts;
    my @todelAlts;
    my @toKeepAlts;
    my $found;

    print("\nDEBUG: vcf_prune\n\n");
    foreach my $variant_to_remove (keys %{$vcf_todelete})
    {
        print("DEBUG: Variant $variant_to_remove\n ");
		if ($tag1==0 and exists($variants1{$variant_to_remove}))
		{
            print("DEBUG: Variant to be removed in vcf1\t\n");
            if ($vcf_todelete->{$variant_to_remove} eq $variants1{$variant_to_remove})
            {
                print("\tDEBUG: removing it in vcf1, $vcf_todelete->{$variant_to_remove} eq $variants1{$variant_to_remove}\n");
                delete $variants1{$variant_to_remove};
            }
            else ## We need to assess if some alleles need to be removed or none of them
            {
                print("\tDEBUG: alleles in vcf1 are not concordant,$vcf_todelete->{$variant_to_remove} vs. $variants1{$variant_to_remove}\n");
                @refAlts=split($OFS,$variants1{$variant_to_remove});
                @todelAlts=split($OFS,$vcf_todelete->{$variant_to_remove});
                if(scalar @refAlts == 2 && scalar @todelAlts == 2) ##Same variant but different alleles
                {
                    print("\t\tDEBUG: different alternative alleles, $vcf_todelete->{$variant_to_remove} vs. $variants1{$variant_to_remove}. Nothing will be removed\n");
                }
                else
                {
                    @toKeepAlts=($refAlts[0]); ##Will contain the set of alternatives that are not in common (will be kept)
                   
                    ##These loops start in 1, since the 0 element is the reference allele
                    for (my $i=1; $i<scalar @refAlts; ++$i)
                    {
                        $found=0;
                        for (my $j=1; $j<scalar @todelAlts; ++$j)
                        {
                            if($refAlts[$i] eq $todelAlts[$j])##This will generate some extra comparisons if pairs have been already found. However, the number of comparisons is very small, I don't think generating extra hashes will save much time
                            {
                                $found=1;
                                last;
                            }
                            
                        }
                        if($found==0)
                        {
                            push(@toKeepAlts,$refAlts[$i]);
                        }
                    }
                    
                    ##If there are alleles that should not be deleted, we update them in the set
                    if(scalar @toKeepAlts >1)
                    {
                        $variants1{$variant_to_remove}=join($OFS,@toKeepAlts);
                        print("\t\tDEBUG: there were a number of alleles that should not be eliminated,$variants1{$variant_to_remove}\n");
                    }
                    else ##There were extra alleles in the list to delete. Thus, we have to eliminate this variant from the problem
                    {
                        print("\t\tDEBUG: there were extra alleles in the reference, but not in the problem. Eliminating variant $variant_to_remove;\n");
                        delete $variants1{$variant_to_remove};
                    }
                }
                
            }
		}
		if ($tag2==0 and exists($variants2{$variant_to_remove}))
		{
            print("DEBUG: Variant to be removed in vcf2\t\n");
            if ($vcf_todelete->{$variant_to_remove} eq $variants2{$variant_to_remove})
            {
                print("\tDEBUG: removing it in vcf2, $vcf_todelete->{$variant_to_remove} eq $variants2{$variant_to_remove}\n");
                delete $variants2{$variant_to_remove};
            }
            else ## We need to assess if some alleles need to be removed or none of them
            {
                print("\tDEBUG: alleles in vcf2 are not concordant, $vcf_todelete->{$variant_to_remove} vs. $variants2{$variant_to_remove}\n");
                @refAlts=split($OFS,$variants2{$variant_to_remove});
                @todelAlts=split($OFS,$vcf_todelete->{$variant_to_remove});
                if(scalar @refAlts == 2 && scalar @todelAlts == 2) ##Same variant but different alleles
                {
                    print("\t\tDEBUG: different alternative alleles, $vcf_todelete->{$variant_to_remove} vs. $variants2{$variant_to_remove}. Nothing will be removed\n");
                }
                else
                { 
                    @toKeepAlts=($refAlts[0]); ##Will contain the set of alternatives that are not in common (will be kept)
                   
                    ##These loops start in 1, since the 0 element is the reference allele
                    for (my $i=1; $i<scalar @refAlts; ++$i)
                    {
                        $found=0;
                        for (my $j=1; $j<scalar @todelAlts; ++$j)
                        {
                            if($refAlts[$i] eq $todelAlts[$j])##This will generate some extra comparisons if pairs have been already found. However, the number of comparisons is very small, I don't think generating extra hashes will save much time
                            {
                                $found=1;
                                last;
                            }
                            
                        }
                        if($found==0)
                        {
                            push(@toKeepAlts,$refAlts[$i]);
                        }
                    }
                    
                    ##If there are alleles that should not be deleted, we update them in the set
                    if(scalar @toKeepAlts >1)
                    {
                        $variants2{$variant_to_remove}=join($OFS,@toKeepAlts);
                        print("\t\tDEBUG: there were a number of alleles that should not be eliminated,$variants2{$variant_to_remove}\n");
                    }
                    else ##There were extra alleles in the list to delete. Thus, we have to eliminate this variant from the problem
                    {
                        print("\t\tDEBUG: there were extra alleles in the reference, but not in the problem. Eliminating variant $variant_to_remove;\n");
                        delete $variants2{$variant_to_remove};
                    }
                }
                
            }
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
# Not genotype specific
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

#Filter out variants in the hash of covN variants. The one staying will be removed from the variants considered in this study
# Not genotype specific
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
    my @refAlts;
    my @todelAlts;
    my @toKeepAlts;
    my $found;

    print("\nDEBUG: vcf_prune_single\n\n");
    foreach my $variant_to_remove (keys %{$vcf_todelete})
    {
#print("DEBUG: Variant $variant_to_remove ");
		if (exists($variants1{$variant_to_remove}))
		{
            print("DEBUG: Variant to be removed in vcf1\n");
            if ($vcf_todelete->{$variant_to_remove} eq $variants1{$variant_to_remove})
            {
                print("\tDEBUG: removing it in vcf1, $vcf_todelete->{$variant_to_remove} eq $variants1{$variant_to_remove}\n");
                delete $variants1{$variant_to_remove};
            }
            else ## We need to assess if some alleles need to be removed or none of them
            {
                print("\tDEBUG: alleles in vcf1 are not concordant, $vcf_todelete->{$variant_to_remove} vs. $variants1{$variant_to_remove}\n");
                @refAlts=split($OFS,$variants1{$variant_to_remove});
                @todelAlts=split($OFS,$vcf_todelete->{$variant_to_remove});
                if(scalar @refAlts == 2 && scalar @todelAlts == 2) ##Same variant but different alleles
                {
                    print("\t\tDEBUG: different alternative alleles, $vcf_todelete->{$variant_to_remove} vs. $variants1{$variant_to_remove}. Nothing will be removed\n");
                }
                else
                {
                    @toKeepAlts=($refAlts[0]); ##Will contain the set of alternatives that are not in common (will be kept)
                   
                    ##These loops start in 1, since the 0 element is the reference allele
                    for (my $i=1; $i<scalar @refAlts; ++$i)
                    {
                        $found=0;
                        for (my $j=1; $j<scalar @todelAlts; ++$j)
                        {
                            if($refAlts[$i] eq $todelAlts[$j])##This will generate some extra comparisons if pairs have been already found. However, the number of comparisons is very small, I don't think generating extra hashes will save much time
                            {
                                $found=1;
                                last;
                            }
                            
                        }
                        if($found==0)
                        {
                            push(@toKeepAlts,$refAlts[$i]);
                        }
                    }
                    
                    ##If there are alleles that should not be deleted, we update them in the set
                    if(scalar @toKeepAlts >1)
                    {
                        $variants1{$variant_to_remove}=join($OFS,@toKeepAlts);
                        print("\t\tDEBUG: there were a number of alleles that should not be eliminated,$variants1{$variant_to_remove}\n");
                    }
                    else ##There were extra alleles in the list to delete. Thus, we have to eliminate this variant from the problem
                    {
                        print("\t\tDEBUG: there were extra alleles in the reference, but not in the problem. Eliminating variant $variant_to_remove;\n");
                        delete $variants1{$variant_to_remove};
                    }
                }
                
            }
		}
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
    my $ref;
    my @aAs;
    my @bAs;
    my %tempAs;

    print("\nDEBUG: vcf_unite_parsed\n\n");
    foreach my $common_variant (keys %common_variants) ##We have to delete possible variants that are not different any more
    {
        if(exists $ref_common_variantsA->{$common_variant})
        {
            if(exists $ref_common_variantsB->{$common_variant})
            {
                #Allele info exists in both. We need to make sure it is the same, or otherwise combine them
                if ($ref_common_variantsA->{$common_variant} eq $ref_common_variantsB->{$common_variant}) #We can get it from either A or B
                {
                    $common_variants{$common_variant}=$ref_common_variantsA->{$common_variant};
                }
                else ##Worst case scenario, we gotta combine them
                {
                    @aAs=split($OFS,$ref_common_variantsA->{$common_variant});
                    @bAs=split($OFS,$ref_common_variantsB->{$common_variant});
                    $ref=shift(@aAs);
                    shift(@bAs);
                    %tempAs=();
                    @tempAs{@aAs,@bAs}=(1) x (scalar @aAs + scalar @bAs);#Union
                    $common_variants{$common_variant}=join($OFS,$ref,keys %tempAs); 
                }
            }
            else #It does not exist in cB, so we are sure we can get it from cA
            {
                $common_variants{$common_variant}=$ref_common_variantsA->{$common_variant};
            }
        }
        else ##If the variant is common but does not exist in cA, we are sure we can get the allele information from cB
        {
            $common_variants{$common_variant}=$ref_common_variantsB->{$common_variant};
        }
        
        if(exists $different_variants{$common_variant})
        {
            print("\tDEBUG: how is this possible?\n");
            ##TODO IMPORTANT UNFINISHED: Need to check when/why this delete can happen (right now I am not able to picture it) and add code for allele consideration in the union of differences (i.e., I can have two different in the same position but different alleles!)
            #delete($different_variants{$common_variant});
        }
    }
    
    foreach my $different_variant (keys %different_variants) ##We have to sync the alleles from the two sides since there is a remote possibility of the same location being private with different alleles in the two sides
    {
        if(exists $ref_different_variantsA->{$different_variant})
        {
            if(exists $ref_different_variantsB->{$different_variant})
            {
                #Allele info exists in both. We need to make sure it is the same, or otherwise combine them
                if ($ref_different_variantsA->{$different_variant} eq $ref_different_variantsB->{$different_variant}) #We can get it from either A or B
                {
                    die "This should definetly not happen\n"
                    #$different_variants{$different_variant}=$ref_differentvariantsA{$different_variant};
                }
                else ##Worst case scenario, we gotta combine them
                {
                    @aAs=split($OFS,$ref_different_variantsA->{$different_variant});
                    @bAs=split($OFS,$ref_different_variantsB->{$different_variant});
                    $ref=shift(@aAs);
                    shift(@bAs);
                    %tempAs=();
                    @tempAs{@aAs,@bAs}=(1) x (scalar @aAs + scalar @bAs);#Union
                    $different_variants{$different_variant}=join($OFS,$ref,keys %tempAs); 
                }
            }
            else #It does not exist in cB, so we are sure we can get it from cA
            {
                $different_variants{$different_variant}=$ref_different_variantsA->{$different_variant};
            }
        }
        else ##If the variant is different but does not exist in cA, we are sure we can get the allele information from cB
        {
            $different_variants{$different_variant}=$ref_different_variantsB->{$different_variant};
        }
       
        ##I will need to make sure that it is not here and if it is work on the alleles 
        if(exists $different_variants{$different_variant})
        {
            print("\tDEBUG: how is this possible?\n");
            ##WORKING HERE: Need to check when/why this delete can happen (right now I am not able to picture it) and add code for allele consideration in the union of differences (i.e., I can have two different in the same position but different alleles!)
            delete($different_variants{$different_variant});
        }
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
    my @aAs;
    my @bAs;
    my @adAs;
    my @bdAs;
    my %tempAs;
    my @cAs;
    my @ctempAs;
    my @candidateAs;
    my $i;
    my $j;
    my $ref;
    my $flag;

    print("\nDEBUG: vcf_intersect_parsed\n\n");
    foreach my $key (keys %{$ref_common_variantsA})
    {
       if(exists $ref_common_variantsB->{$key})
       {
           #Allele info exists in both. We need to make sure it is the same, or otherwise choose the commons. If no commons, is not in the intersection
           if ($ref_common_variantsA->{$key} eq $ref_common_variantsB->{$key}) #We can get it from either A or B
           {
               $common_variants{$key}=$ref_common_variantsA->{$key};
           }
           else ##Only the commons, if there are commons
           {
                @aAs=split($OFS,$ref_common_variantsA->{$key});
                @bAs=split($OFS,$ref_common_variantsB->{$key});
                @cAs=(shift(@aAs)); #We initialize the common array with the reference
                shift(@bAs);
                for ($i=0; $i<scalar @aAs; ++$i)
                {
                    for ($j=0; $j<scalar @bAs; ++$j)
                    {
                        if($aAs[$i] eq $bAs[$j])
                        {
                            push(@cAs,$aAs[$i]);
                            last;
                        }
                    }
                }
                if (scalar @cAs > 1)
                {
                    $common_variants{$key}=join($OFS,@cAs); 
                }
           }
       }
    }
    
    my %different_variants=(%{$ref_different_variantsA},%{$ref_different_variantsB},%{$ref_common_variantsA},%{$ref_common_variantsB});#Union of all variants

    foreach my $variant (keys %different_variants)
    {
        @aAs=();
        @bAs=();
        @adAs=();
        @bdAs=();

        if(exists $ref_different_variantsA->{$variant})
        {
            @adAs=split($OFS,$ref_different_variantsA->{$variant});
            $ref=shift(@adAs);
        }
        if(exists $ref_different_variantsB->{$variant})
        {
            @bdAs=split($OFS,$ref_different_variantsB->{$variant});
            $ref=shift(@bdAs);
        }
        if(exists $ref_common_variantsA->{$variant})
        {
            @aAs=split($OFS,$ref_common_variantsA->{$variant});
            $ref=shift(@aAs);
        }
        if(exists $ref_common_variantsB->{$variant})
        {
            @bAs=split($OFS,$ref_common_variantsB->{$variant});
            $ref=shift(@bAs);
        }

        %tempAs=();
        @tempAs{@aAs,@bAs,@adAs,@bdAs}=(1) x (scalar @aAs + scalar @bAs + scalar @adAs + scalar @bdAs);#Union of all alternate alleles

        @candidateAs=keys %tempAs;

        if (exists $common_variants{$variant}) #if it is in common, we can only accept the alternative alleles that are not in common, or delete this variant
        {
            @cAs=($ref);
            @ctempAs=$common_variants{$variant};
            shift(@ctempAs);
            for ($i=0; $i<scalar @candidateAs; ++$i)
            {
                $flag=0;
                for($j=0; $j<scalar @ctempAs; ++$j)
                {
                    if($candidateAs[$i] eq $ctempAs[$j])
                    {
                        $flag=1;
                        last;
                    }
                }
                if($flag==0)
                {
                    push(@cAs,$candidateAs[$i]);
                }
            }
            if(scalar @cAs >1) ##Despite having some common alleles for this position, we also have some different
            {
                $different_variants{$variant}=join($OFS,@cAs);
            }
            else ##All are common. I need to delete this variant
            {
                delete $different_variants{$variant};
            }
        }
        else #if it is not in common, this is a good different variant
        {
            $different_variants{$variant}=join($OFS,$ref,@candidateAs);
        }
    }

    my $n_selvariants=scalar keys %common_variants;
    my $n_filtvariants=scalar keys %different_variants;
    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
    return (\%common_variants,\%different_variants);
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
    my $value;
	for (my $i=1; $i<scalar @{$array};$i++) ###The first line is the header
	{
        chomp(${$array}[$i]);
		$key=$value=${$array}[$i];
		$key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
        $value=~s/^[^\t]+\t[^\t]+\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2/;
		$hash{$key}=$value;
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
# This subrutine is genotype compatible: it prints only the information relative to
# the genotypes present in the hash, even if several where present in the VCF file
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
    my @newline;
    my @helper; ##I know I know, I am microoptimizing and making this unreadable
    my @helper2;
    my $helpstring;
    my @contline;
    my @vcfAs;
    my @hashAs;
    my $vcfValue;
    my @validAs;
    my $i;
    my $j;
    my $k;
     
    #Copying the header and adding a new line with filtering info
    #Then adding the variants that are present in the hash.
    for($i=0;$i< scalar @icontent; ++$i)
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
                if($icontent[$i]=~m/^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+,[^\t]+\t.*/) #If there are several alternative alleles, we need to do extra work
                {
                    chomp($icontent[$i]);
                    @contline=split("\t", $icontent[$i]);
                    @vcfAs=split(",",$contline[$4]);
                    $vcfValue=join($OFS, $contline[$3],@vcfAs);
                    if ($vcfValue eq $hash{$key}) ##All alleles included, just print the line
                    {
                        print($OFILE $icontent[$i],"\n");
                    }
                    else #We need to detect which alleles need to be printed
                    {
                        @validAs=();
                        @hashAs=split($OFS,$hash{$key});
                        shift(@hashAs);
                        for($j=0; $j<scalar @vcfAs; ++$j)
                        {
                            for($k=0; $k<scalar @hashAs; ++$k)
                            {
                                if($vcfAs[$j] eq $hashAs[$k])
                                {
                                    push(@validAs,$j);
                                    last;
                                }
                            }
                        }
                        if(scalar @vcfAs == scalar @hashAs) ##Valid alleles in the VCF, need to print
                        {
                            @newline=();
                            push(@newline,@contline[0 .. 3]); #CHROM  POS ID  REF
                            push(@newline,join(",",@vcfAs[@validAs])); #ALT
                            push(@newline,@contline[5,6]); #QUAL    FILTER
                            @helper=split(";",$contline[7]); #INFO, format: X..X=???,???...???;
                            for ($j=0; $j< scalar @helper; ++$j) #We traverse the INFO fields, ;-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                            {
                                if($helper[$j] =~ m/,/)
                                {
                                    $helpstring=$helper[$j];
                                    $helpstring=~s/^(.+=).+$/$1/; #header of this field
                                    $helper[$j]=~s/^.+=(.+)$/$1/; #contents
                                    @helper2=split(",",$helper[$j]);
                                    $helper[$j]=join(",",@helper2[@validAs]);
                                    $helper[$j]=$helpstring.$helper[$j];
                                }
                            }
                            push(@newline,join(";",@helper));## Add reconstructed info field back
                            push(@newline,$contline[8]); #FORMAT
                            
                            #SAMPLES. This should be just one in this approach, but I am generalizing it
                            for ($j=9; $j<scalar @contline; ++$j)
                            {
                                @helper=split(":",$contline[$j]); #SAMPLE, format: fields :-separated, some, allele separated (,).
                                for ($k=0; $k< scalar @helper; ++$k) #We traverse the SAMPLE fields, :-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                                {
                                    if($helper[$k] =~ m/,/)
                                    {
                                        @helper2=split(",",$helper[$k]);
                                        $helper[$k]=join(",",@helper2[@validAs]);
                                    }
                                }
                                push(@newline,join(":",@helper));## Add reconstructed SAMPLE field back
                                 
                            }
                            print($OFILE join("\t",@newline),"\n");
                        }
                        else
                        {
                            die "ERROR in write_variant_vcf. The alleles present in the hash are not present in the VCF file\n\n";
                        } #Else we don't do anything, this hash variant will be removed but not printed
                    }
                }
                else #Otherwise we are good to go
                {
                    print($OFILE $icontent[$i],"\n");
                }
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
# present in the two VCF files. The header comes also from the first VCF file. ONLY one sample information, comming from the first VCF file with higher priority.
# This subrutine is compatible with genotypes. If a given genotype is present in the two
# VCF files, the information will come from the first VCF file. If different genotypes
# come from different VCF files, all information comes from the first file except the one 
# relative to the alleles that are not present in that file. THIS IS SOME KIND OF FRANKENSTEIN'S
# MONSTER VCF FILE. IT WILL DEFINETLY CONTAIN INCONSISTENCIES. USE IT ONLY AS A LIST OF VARIANTS
# ##################################################################################
sub write_variant_2vcf
{
    my ($ref_hash,$filename,$vcf,$vcf2,$comment,$filter)=@_;
    my $OFILE=openwrite_vcf_refname("$filename",$filter);
    my %outcontent;
    (! defined $OFILE) and return;
    my %vcf2_hash=parse_vcf_complete($vcf2);
    open(my $IFILE, $vcf);
    open(my $IFILE2, $vcf2);
    my @icontent= <$IFILE>;
    my @icontent2= <$IFILE2>;
    close($IFILE);
    close($IFILE2);
    my $flag=0;
    my %hash=%{$ref_hash};
    my $key;
    my $value;
    my @contline;
    my @contline_vcf2;
    my @ids_contline;
    my @newline;
    my @helper;
    my @helper2;
    my @helper21;
    my @helper22;
    my @helper3;
    my $helpstring;
    my $helpstring2;
    my @vcfAs;
    my @vcfAs2;
    my @validAs;
    my @hashAs;
    my $vcfValue;
    my $j;
    my $k;
    my $l;
     
    print("\nDEBUG: write_variant_2vcf\n\n");

    ##Description of the main strategy (this subrutine is a mess):
    #   for each line of vcf1
    #       if the vcf matches the alleles in the hash, print
    #       elsif, if they are a superset of the alleles in hash, print the needed
    #       elsif they are not enough, get from vcf2
    #       else , do nothing, this will get parsed from vcf2


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
        else ##Most lines
        {
            $key=$icontent[$i];
            chomp($key);
            $value=$key;
            $key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
            $value=~s/^[^\t]+\t[^\t]+\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2/;
            #print("DEBUG: Key $key\n");
            if(exists $hash{$key})
            {
                if($value eq $hash{$key})
                {
                    ##This file contains all needed information
                    $outcontent{$key}=$icontent[$i];
                    delete($hash{$key});
                    delete($vcf2_hash{$key});
                }
                else
                {
                    chomp($icontent[$i]);
                    @contline=split("\t", $icontent[$i]);
                    @vcfAs=split(",",$contline[$4]);
                    @validAs=();
                    @hashAs=split($OFS,$hash{$key});
                    shift(@hashAs);
                    for($j=0; $j<scalar @vcfAs; ++$j)
                    {
                        for($k=0; $k<scalar @hashAs; ++$k)
                        {
                            if($vcfAs[$j] eq $hashAs[$k])
                            {
                                push(@validAs,$j);
                                last;
                            }
                        }
                    }
                    if (scalar @validAs == scalar @hashAs) ##All needed alleles in vcf1, no need to use vcf2
                    {
                        @newline=();
                        push(@newline,@contline[0 .. 3]); #CHROM  POS ID  REF
                        push(@newline,join(",",@vcfAs[@validAs])); #ALT
                        push(@newline,@contline[5,6]); #QUAL    FILTER
                        @helper=split(";",$contline[7]); #INFO, format: X..X=???,???...???;
                        for ($j=0; $j< scalar @helper; ++$j) #We traverse the INFO fields, ;-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                        {
                            if($helper[$j] =~ m/,/)
                            {
                                $helpstring=$helper[$j];
                                $helpstring=~s/^(.+=).+$/$1/; #header of this field
                                $helper[$j]=~s/^.+=(.+)$/$1/; #contents
                                @helper2=split(",",$helper[$j]);
                                $helper[$j]=join(",",@helper2[@validAs]);
                                $helper[$j]=$helpstring.$helper[$j];
                            }
                        }
                        push(@newline,join(";",@helper));## Add reconstructed info field back
                        push(@newline,$contline[8]); #FORMAT
                        
                        #SAMPLES. This should be just one in this approach, but I am generalizing it
                        for ($j=9; $j<scalar @contline; ++$j)
                        {
                            @helper=split(":",$contline[$j]); #SAMPLE, format: fields :-separated, some, allele separated (,).
                            for ($k=0; $k< scalar @helper; ++$k) #We traverse the SAMPLE fields, :-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                            {
                                if($helper[$k] =~ m/,/)
                                {
                                    @helper2=split(",",$helper[$k]);
                                    $helper[$k]=join(",",@helper2[@validAs]);
                                }
                            }
                            push(@newline,join(":",@helper));## Add reconstructed SAMPLE field back
                             
                        }
                        $outcontent{$key}=join("\t",@newline);
                        delete($hash{$key});
                        delete($vcf2_hash{$key});
                    }
                    elsif (scalar @validAs != 0) ##Combine alleles from vcf1 and vcf2
                    {
                        @contline_vcf2=split("\t",$vcf2_hash{$key});
                        scalar @contline_vcf2 != scalar @contline and die "Error, the number of VCF columns is not the same in the two VCF files\n\n";
                        @vcfAs2=split(",",$contline_vcf2[4]);
                        
                        for($j=0; $j<scalar @hashAs; ++$j)
                        {
                            $helpstring=0;
                            for($k=0; $k<scalar @vcfAs; ++$k)
                            {
                                if($hashAs[$j] eq $vcfAs[$k])
                                {
                                    $ids_contline[$j]=1;
                                    push(@validAs,$k);
                                    $helpstring=1;
                                    last;
                                }
                            }
                            if($helpstring==0)
                            {
                                for($k=0; $k<scalar @vcfAs2; ++$k)
                                {
                                    if($hashAs[$j] eq $vcfAs2[$k])
                                    {
                                        $ids_contline[$j]=2;
                                        push(@validAs,$k);
                                        $helpstring=1;
                                        last;
                                    }
                                }
                            }
                            if($helpstring==0)
                            {
                                die "Error in write_variant_2vcf";
                            }
                        }
                        
                        @newline=();
                        push(@newline,@contline[0 .. 3]); #CHROM  POS ID  REF
                        @helper=();
                        for ($j=0; $j< scalar @validAs; ++$j)
                        {
                            push(@helper,(split(",",$ids_contline[$j]==1?$contline[4]:$contline_vcf2[4]))[$validAs[$j]]);
                        }
                        push(@newline,join(",",@helper)); #ALT
                        push(@newline,@contline[5,6]); #QUAL    FILTER
                        @helper=split(";", $contline[7]); #INFO
                        @helper2=split(";", $contline_vcf2[7]); #INFO_VCF2
                        scalar @helper != scalar @helper2 and die "Error in write_variant_2vcf. The number of elements of the INFO field is not the same in the two VCF files\n\n";
                        
                        for ($j=0; $j< scalar @helper; ++$j) #We traverse the INFO fields, ;-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                        {
                            if($helper[$j] =~ m/,/ || $helper2[$j] =~ m/,/)
                            {
                                $helpstring=$helper[$j];
                                $helpstring2=$helper2[$j];
                                
                                $helpstring=~s/^(.+=).+$/$1/; #header of this field
                                $helpstring2=~s/^(.+=).+$/$1/; #header of this field
                                $helpstring ne $helpstring2 and die "Error in write_variant_2vcf. The INFO fields of vcf1 and 2 are not the same/in the same order";
                                $helper[$j]=~s/^.+=(.+)$/$1/; #contents
                                $helper2[$j]=~s/^.+=(.+)$/$1/; #contents

                                @helper21=split(",",$helper[$j]);
                                @helper22=split(",",$helper2[$j]);
                                @helper3=();
                                for ($k=0; $k< scalar @validAs; ++$k)
                                {
                                    push(@helper3,$ids_contline[$j]==1?$helper21[$validAs[$k]]:$helper22[$validAs[$k]]);
                                }
                                $helper[$j]=$helpstring.join(",",@helper3);
                            }
                        }
                        push(@newline,join(";",@helper));## Add reconstructed info field back
                        push(@newline,$contline[8]); #FORMAT
                        
                        #SAMPLES. This should be just one in this approach, but I am generalizing it
                        for ($j=9; $j<scalar @contline; ++$j)
                        {
                            @helper=split(":",$contline[$j]); #SAMPLE, format: fields :-separated, some, allele separated (,).
                            @helper2=split(":",$contline_vcf2[$j]); #SAMPLE, format: fields :-separated, some, allele separated (,).
                            for ($k=0; $k< scalar @helper; ++$k) #We traverse the SAMPLE fields, :-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                            {
                                if($helper[$k] =~ m/,/ || $helper2[$k] =~ m/,/)
                                {
                                    @helper21=split(",",$helper[$k]);
                                    @helper22=split(",",$helper2[$k]);
                                    @helper3=();
                                    for ($l=0; $l<scalar @validAs; ++$l)
                                    {
                                        push(@helper3,$ids_contline[$j]==1?$helper21[$validAs[$l]]:$helper22[$validAs[$l]]);
                                    }
                                    $helper[$k]=join(",",@helper3);
                                }
                            }
                            push(@newline,join(":",@helper));## Add reconstructed SAMPLE field back
                             
                        } 
                        $outcontent{$key}=join("\t",@newline);
                        delete($hash{$key});
                        delete($vcf2_hash{$key});
                    }
#                    else ##Don't do anyting, this variable will be dealt with after the loop, using vcf2
#                    {
#                    }
                }
            }
            if(scalar keys %hash == 0)
            {
                last;
            }
            
        }

    }
    if (scalar keys %hash !=0) #Not present at all in vcf1 (the needed alleles), needs to be printed from vcf2
    {
        foreach $key (keys %hash)
	    {
            exists $vcf2_hash{$key} or die "Error in write_variant_2vcf. Variant that was not cleared using VCF1 is not present in VCF2\n\n";
    		
			$value=$vcf2_hash{$key};
            $value=~s/^[^\t]+\t[^\t]+\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2/;
       		#print("DEBUG: Key $key\n");
            if($value eq $hash{$key}) ##This file contains all needed information as needed
            {
                $outcontent{$key}=$vcf2_hash{$key};
            }
            elsif($vcf2_hash{$key}=~m/^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+,[^\t]+\t.*/) ##This file contains all needed information, but some alleles need to be discarded
            {
                @contline=split("\t",$vcf2_hash{$key});
                @vcfAs=split(",",$contline[$4]);
                @validAs=();
                @hashAs=split($OFS,$hash{$key});
                shift(@hashAs);
                for($j=0; $j<scalar @vcfAs; ++$j)
                {
                    for($k=0; $k<scalar @hashAs; ++$k)
                    {
                        if($vcfAs[$j] eq $hashAs[$k])
                        {
                            push(@validAs,$j);
                            last;
                        }
                    }
                }
                scalar @validAs != scalar @hashAs and die "ERROR: Variant only present in VCF2 does not have all needed alleles present in VCF2. Where did these alleles come from?\n";
                @newline=();
                push(@newline,@contline[0 .. 3]); #CHROM  POS ID  REF
                push(@newline,join(",",@vcfAs[@validAs])); #ALT
                push(@newline,@contline[5,6]); #QUAL    FILTER
                @helper=split(";",$contline[7]); #INFO, format: X..X=???,???...???;
                for ($j=0; $j< scalar @helper; ++$j) #We traverse the INFO fields, ;-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                {
                    if($helper[$j] =~ m/,/)
                    {
                        $helpstring=$helper[$j];
                        $helpstring=~s/^(.+=).+$/$1/; #header of this field
                        $helper[$j]=~s/^.+=(.+)$/$1/; #contents
                        @helper2=split(",",$helper[$j]);
                        $helper[$j]=join(",",@helper2[@validAs]);
                        $helper[$j]=$helpstring.$helper[$j];
                    }
                }
                push(@newline,join(";",@helper));## Add reconstructed info field back
                push(@newline,$contline[8]); #FORMAT
                
                #SAMPLES. This should be just one in this approach, but I am generalizing it
                for ($j=9; $j<scalar @contline; ++$j)
                {
                    @helper=split(":",$contline[$j]); #SAMPLE, format: fields :-separated, some, allele separated (,).
                    for ($k=0; $k< scalar @helper; ++$k) #We traverse the SAMPLE fields, :-separated. The ones that have several fields (,-separated) are broken up and only the valid indexes are chosen.
                    {
                        if($helper[$k] =~ m/,/)
                        {
                            @helper2=split(",",$helper[$k]);
                            $helper[$k]=join(",",@helper2[@validAs]);
                        }
                    }
                    push(@newline,join(":",@helper));## Add reconstructed SAMPLE field back
                     
                }
                $outcontent{$key}=join("\t",@newline);
            }
            else
            {
                die "ERROR: in write_variant_2vcf. Variant not present in VCF1 and incomplete in VCF2. Where did the other alleles come from?";
            }
        }
    }
    
    ##Sort outcontent and print to OFILE

    my @int_v;

    foreach $key (nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %outcontent)
    {
        print($OFILE $outcontent{$key},"\n");
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
    ##NAB constants
    foreach my $cond (@NABfiltering_conditions)
    {
        $constExeCondContentRefs{"${exe_condition}_${cond}"}=vcf_covN_filter($constExeCondContentRefs{$covNname},$cond);
    }

    $constExeCondContentRefs{$Aexecondname}=parse_vcf($Aexecondname);
    $constExeCondContentRefs{$Bexecondname}=parse_vcf($Bexecondname);
    
    ##PAF constants
    $constExeCondContentRefs{"${exe_condition}_PAF"}=getPAFdata($constExeCondContentRefs{$Aexecondname},$constExeCondContentRefs{$Bexecondname});
  
    ##TODO: Do I want to parallelize this independently? Same question about the benefits and costs of forking
    foreach my $cond (@popAFfiltering_conditions)
    {
        $constExeCondContentRefs{"${exe_condition}_${cond}"}=filt_PAF($constExeCondContentRefs{"${exe_condition}_PAF"},$cond);
    } 
    
    $nofiltResultsRef{$exe_condition}=[scalar keys %{$constExeCondContentRefs{$AcovBname}}, scalar keys %{$constExeCondContentRefs{$BcovBname}}, scalar keys %N, scalar keys %{vcf_prune_single($constExeCondContentRefs{$Aexecondname},\%N)}, scalar keys %{vcf_prune_single($constExeCondContentRefs{$Bexecondname},\%N)}];
}

##Returns a hash ref to the population allele frequency information for all variants present in the input reference hashes 
sub getPAFdata
{
    my $tabix = Bio::DB::HTS::Tabix->new( filename =>$GNOMAD );
    my %outdata; #key: CHROM${OFS}POS${OFS}ALT value: [$tref,$tfilt,$taf];
    my @keyOutPAFs;
    my @keys;
    my $ref;
    my $parsedPAFS;
    my %obtainedPAFs; #key = CHROM${OFS}POS #value = {ALT => [$tref, $tfilt, $taf]}
    my %neededKeys;
    my $tabix_iter;
    my $chr;
    my $nstart;
    my $line;
    my $flag;
    my $i;
    my ($tstart,$talt,$tref,$tfilt,$taf);
    my @alts;

    foreach my $ref_hash (@_)
    {
        foreach my $key (keys %{$ref_hash}) ##Foreach key
        {
            #key = CHROM${OFS}POS
            #$ref_hash->{$key} = REF${OFS}ALT1${OFS}ALT2...${OFS}ALTN
            unless (exists $obtainedPAFs{$key})
            {
                ##Obtain gnomAD data for this genomic position if it has not been obtained before. MULI-SNVs are represented in different lines
                ($chr,$nstart)=split($OFS,$key);
                $tabix_iter=$tabix->query("$chr:$nstart-".($nstart+1));
                $flag=0;
                if(defined $tabix_iter)
                {
                    while($line=$tabix_iter->next)
                    {
                        #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  S1 ... SN
                        $line =~ s/^[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t[^\t]+AF=([^;]+).*$/$1\t$3\t$2\t$4\t$5/;
                        ($tstart,$talt,$tref,$tfilt,$taf)=split("\t",$line);
                        if($flag==0) ##Need to generate the new hash
                        {
                            my %thisposdata;
                            $thisposdata{$talt}=[$tref,$tfilt,$taf];
                            $obtainedPAFs{"$chr$OFS$nstart"}=\%thisposdata;
                            $flag=1;
                        }
                        else
                        {
                            $obtainedPAFs{"$chr$OFS$nstart"}->{$talt}=[$tref,$tfilt,$taf];
                        }
                    }
    
                }
                else
                {
                    warn "There is no data for $chr:$nstart-".($nstart+1);
                }
            }
            @alts=split($OFS,$ref_hash->{$key});
            $ref=shift(@keys);
            @keyOutPAFs= map {$key.$OFS.$_} @alts;
            
            foreach ($i=0; $i<scalar @keyOutPAFs; ++$i)
            {
                unless (exists $outdata{$keyOutPAFs[$i]}) ##If it exists we don't need to do anything
                {
                    if(exists $obtainedPAFs{$key} && exists $obtainedPAFs{$key}->{$alts[$i]})
                    {
                        $outdata{$keyOutPAFs[$i]}=$obtainedPAFs{$key}->{$alts[$i]};
                    }
                    else
                    {
                        $outdata{$keyOutPAFs[$i]}=["NA","NA","NA"];
                    }
                }
            }
        }##foreach key
    }##foreach hash
    return \%outdata;
}


##This subrutine gets an array of PAF conditions and generates the frequency of private A, B and common C variants with PAF<=condition for each condition
sub getPAFStats
{
    my ($ref_conditions, $ref_pAFfiltN, $ref_A, $ref_B, $ref_C)=@_;
    my @output;
    my @filt_values= map {$_=~s/^--max_pAF://} @{$ref_conditions};
    my $ngood;
    my $ntotal;
    my $pos_key;
    my $comp_key;
    my @alts;
    #ref_pAFfiltN; #key: CHROM${OFS}POS${OFS}ALT value: [$tref,$tfilt,$taf];

    foreach my $value (@filt_values)
    {
        #PrivA
        $ngood=0;
        $ntotal=0;
        foreach $pos_key (keys %{$ref_A})
        {
            @alts=split($OFS,$ref_A->{$pos_key});
            shift(@alts);
            foreach $comp_key (map {$pos_key.$OFS.$_} @alts)
            {
                if(exists $ref_pAFfiltN->{$comp_key})
                {
                    if($ref_pAFfiltN->{$comp_key}->[2] <= $value)
                    {
                        ++$ngood;
                    }
                }
                else
                {
                    ++$ngood;
                    warn "No population data for $comp_key\n";
                }
                ++$ntotal;
            }
        }
        push(@output,$ngood*100.0/$ntotal);

        #PrivB
        $ngood=0;
        $ntotal=0;
        foreach $pos_key (keys %{$ref_B})
        {
            @alts=split($OFS,$ref_B->{$pos_key});
            shift(@alts);
            foreach $comp_key (map {$pos_key.$OFS.$_} @alts)
            {
                if(exists $ref_pAFfiltN->{$comp_key})
                {
                    if($ref_pAFfiltN->{$comp_key}->[2] <= $value)
                    {
                        ++$ngood;
                    }
                }
                else
                {
                    ++$ngood;
                    warn "No population data for $comp_key\n";
                }
                ++$ntotal;
            }
        }
        push(@output,$ngood*100.0/$ntotal);

        #Common
        $ngood=0;
        $ntotal=0;
        foreach $pos_key (keys %{$ref_C})
        {
            @alts=split($OFS,$ref_C->{$pos_key});
            shift(@alts);
            foreach $comp_key (map {$pos_key.$OFS.$_} @alts)
            {
                if(exists $ref_pAFfiltN->{$comp_key})
                {
                    if($ref_pAFfiltN->{$comp_key}->[2] <= $value)
                    {
                        ++$ngood;
                    }
                }
                else
                {
                    ++$ngood;
                    warn "No population data for $comp_key\n";
                }
                ++$ntotal;
            }
        }
        push(@output,$ngood*100.0/$ntotal);
    } ##Foreach filter

    return \@output;
}

##Subrutine to generate the header of the PAF stats for the output file
##Inputs: two array references, to a) names of the different output columns
## b) array of PAF filtering conditions
## Return: array of column headers that combines the names with the filtering values
sub makeHeaderPAF
{
    my ($ref_names,$ref_conditions)=@_;
    my @outcolumns;
    foreach my $condition (@{$ref_conditions})
    {
        $condition =~s/--[^:]+://;
        foreach my $name (@{$ref_names})
        {
            push(@outcolumns,$name.$condition);
        }
    }
    return @outcolumns;
}

##Subrutine that generates a list of variants that are valid (<= condition) considering their population allele frequency
sub filt_PAF
{
    my ($ref_data,$condition)=@_;
    $condition =~s/--[^:]+://; ##Get the maximum population allele frequency
    my %outdata;
    foreach my $key (keys %{$ref_data})
    {
        if($ref_data->{$key}[2]<=$condition)
        {
            $outdata{$key}=$ref_data->{$key};
        }
    }

    return \%outdata;
}

##Filters two hashes of variants using a hash of variants that are valid (this will be a superset of the variants that will be valid)  and generates stats like vcf_prune
sub filter_with_PAF
{
    my ($ref_common, $ref_different, $ref_PAF, $ref_stats)=@_;
    my (%out_common, %out_different);
    my $ref;
    my $ckey;
    my $alt;
    my @alts;
    my @validAlts;
    
    ##The number of ref_PAF variants will be much larger than ref_common or ref_different. Thus, I will loop along those two, instead of the first.
    foreach my $variant (keys %{$ref_common})
    {
        @alts=split($OFS,$ref_common->{$variant});
        $ref=shift(@alts);
        @validAlts=($ref);
        foreach $alt (@alts)
        {
            $ckey=$variant.$OFS.$alt;
            if(exists $ref_PAF->{$ckey}) ##I don't need to check the value of this since this list of variants has been pre-made
            {
                push(@validAlts,$alt);
            }
        }
        if(scalar @validAlts > 1)
        {
            $out_common{$variant}=join($OFS,@validAlts);
        }
    }
    foreach my $variant (keys %{$ref_different})
    {
        @alts=split($OFS,$ref_different->{$variant});
        $ref=shift(@alts);
        @validAlts=($ref);
        foreach $alt (@alts)
        {
            $ckey=$variant.$OFS.$alt;
            if(exists $ref_PAF->{$ckey}) ##I don't need to check the value of this since this list of variants has been pre-made
            {
                push(@validAlts,$alt);
            }
        }
        if(scalar @validAlts > 1) ##The 1 is the reference, we need at least one valid alternative 
        {
            $out_different{$variant}=join($OFS,@validAlts);
        }
    }
    	
    my $n_selvariants=scalar keys %out_common;
	my $n_filtvariants=scalar keys %out_different;
    if($n_selvariants+$n_filtvariants==0)
    {
         @{ $ref_stats }=(0,$n_selvariants,$n_selvariants+$n_filtvariants);
    }
    else
    {   	
	    @{ $ref_stats }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants,$n_selvariants+$n_filtvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
	}
    return (\%out_common, \%out_different);
}
