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
our $private="-p private";
our $sep_param="#";
our $sep_value=":";
our $OFS=",";
our $FS=",";
our $variant_caller="platypus";
our $variant_calling_sh;
our $helper_sh="HeterAnalyzer_multiple_helper.sh";
our $helper_pl="HeterAnalyzer_multiple.pl";
our $tstv_sh="tstv.sbatch";
#our $annotation_sh="annovar.sh";
our $n_cores=1;
our $qsub="sbatch ${private} -N 1 -c ";
our $qsub_noparallel="sbatch ${private} -N 1 -n 1 -c 1";
our $qstat="qstat";
our $sed='sed "s/Submitted batch job \(.*\)/\1/"';
our $sleep=60;
######################################################

##IO Variables
######################################################
my $execond_inputfile="";
my $filtercond_inputfile="";
my $NABfiltercond_inputfile1="",
my $NABfiltercond_inputfile2="",
my $output_dir="vcf_outputdir";
my $input_file="";
#my $normal_bam="";
#my $sample1_bam="";
#my $sample2_bam="";
my $SCRIPTSVCF_DIR=$ENV{'SCRIPTSVCF_DIR'};

#Flags
my $help;
my $usage="\nUsage: $0 [options] -i inputfile\n\n\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the NAB sample\n\t--NABfilt_cond_inputfile2 : input file for the secondary filtering options of the NAB sample (OR filter implemented in a dirty way)\n\t--covaltB_cond_inputfile : input file for the filtering taking into account characteristics of the unfiltered in the comparison\n\t--output_dir : output directory for vcf files\n\t--n_cores : number of cores to execute some steps in parallel\n\t\n\n";
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
    'output_dir=s' => \$output_dir,
    'input_file|i=s' => \$input_file,
    'n_cores=i' => \$n_cores,
    'help|h' => \$help,
                )) or (($input_file eq "") || $help) and die $usage;

$qsub.=$n_cores;

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

my $oefile=Cwd::abs_path($execond_inputfile);
my $offile=Cwd::abs_path($filtercond_inputfile);
my $onfile=Cwd::abs_path($NABfiltercond_inputfile1);
my $onfile2=Cwd::abs_path($NABfiltercond_inputfile2);
my $ifile=Cwd::abs_path($input_file);
my $ocfile=Cwd::abs_path($covBfiltercond_inputfile);

mkdir $output_dir;
chdir $output_dir or die "The output directory $output_dir is not accesible";
$output_dir=Cwd::abs_path($output_dir);

if (-f "$SCRIPTSVCF_DIR/$variant_caller.sh")
{
    $variant_calling_sh="$SCRIPTSVCF_DIR/$variant_caller.sh";
}
else
{
    die "Error, the sh file for the variant caller $variant_caller can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR is indicating the folder with this package of scripts\n"
}

if (-f "$SCRIPTSVCF_DIR/$helper_sh" && -f "$SCRIPTSVCF_DIR/$helper_pl")
{
    $helper_sh="$SCRIPTSVCF_DIR/$helper_sh";
    $helper_pl="$SCRIPTSVCF_DIR/$helper_pl";
}
else
{
    die "Error, the files for the secondary analysis of the data, $helper_sh and $helper_pl can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR is indicating the folder with this package of scripts\n"
}

if (-f "$SCRIPTSVCF_DIR/$tstv_sh")
{
    $tstv_sh="$SCRIPTSVCF_DIR/$tstv_sh";
}
else
{
  die "Error, the sh file to get the TsTv data, $tstv_sh can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR is indicating the folder with this package of scripts\n"
}

my $vcf_filt_exe="$SCRIPTSVCF_DIR/vcf_filtering.pl";

if (! -f "$vcf_filt_exe")
{
    die "The executable vcf_filtering.pl can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR indicates the directory with this package of scripts\n";
}

my $postanalyzer_exe="$SCRIPTSVCF_DIR/postHeterAnalyzer_multiple.sh";
if (! -f $postanalyzer_exe)
{
    die "The executable postanalyzer.exe can't be located. Please, make sure that the environment variable SCRIPTSVCF_DIR indicates the directory with this package of scripts\n";
}

##Input file parsing
####################

my %cases;
my %datafiles;
my $normalfile="";
open(my $INPUT,$input_file) or die "Error opening input file $input_file\n$usage";
my @acases=<$INPUT>;
close($INPUT);
my $name;
my $normal;
my $a;
my $b;

foreach my $case (@acases)
{
    ($name,$normal, $a, $b) = split (" ",$case);
    $cases{$name}=[$normal,$a,$b];

##EDIT Why did I do this? I don't think this is right :S
#    if ($normalfile eq "")
#    {
#        $normalfile= $normal;
#    }
#    ($normal ne $normalfile) and die "Error, all normal files must be the same\n";

    $datafiles{$a}=1;
    $datafiles{$b}=1;
}

## Main conditions loop
#######################

my $name_condition;
my @exe_conditions;
my $bamfiles;
my @filtering_conditions;
my @NABfiltering_conditions1; 
my @NABfiltering_conditions2;
my @NABfiltering_conditions;
my @covBfiltering_conditions;

##Generate all combinations of proposed values for execution and filtering parameters
combs(0,"",\@exe_parameters,\@exe_param_values,\@exe_conditions);
combs(0,"",\@filtering_parameters,\@filtering_param_values,\@filtering_conditions);
combs(0,"",\@NABfiltering_parameters1,\@NABfiltering_param_values1,\@NABfiltering_conditions1);
combs(0,"",\@NABfiltering_parameters2,\@NABfiltering_param_values2,\@NABfiltering_conditions2);
combs(0,"",\@covBfiltering_parameters,\@covBfiltering_param_values,\@covBfiltering_conditions);

our %job_ids;

mkdir("e_logs");
mkdir("o_logs");

### Default calling parameters. This may be changed in the future
#########################################################

my $job_id;
my $normalvcf=basename($normalfile);
my $vcfname=$normalvcf;
$vcfname=~s/.bam//;
$normalvcf=~s/.bam/.vcf/;
print("Unfiltered variant calling detection/execution:\n");	
if(! -f $normalvcf)
{
	$bamfiles=$normalfile;
	my $exe_condition="";
    $job_id=submit_job_name($vcfname,"$qsub $variant_calling_sh $bamfiles $vcfname.vcf ${vcfname}_platypus.log $exe_condition");
    print("\tNormal tissue variant calling submited with job_id $job_id\n");
}
else
{
	print("\tNormal tissue variant calling already present, skipping it\n");
}

my @afiles=sort keys %datafiles; #If we do not sort them NAB files will not have always the same order between runs
my $fjob_id="";
my $actual_exe_conditions;

for (my $i=0; $i<scalar(@afiles); ++$i)
{
    $bamfiles=$afiles[$i];
    $vcfname=basename($bamfiles);
    $vcfname=~s/.bam//;
    unless (-f "$vcfname.vcf") ##Variant calling for each cancer file
    {
        $job_id=submit_job_name($vcfname,"$qsub $variant_calling_sh $bamfiles $vcfname.vcf ${vcfname}_platypus.log"); ##Defaults do not generate dependencies
        print("\tSample $vcfname variant calling submited with job_id $job_id\n");
    }
    else
    {
        print("\tSample $vcfname variant calling already present. Skipping\n");
    }

    foreach my $exe_condition (@exe_conditions)
    {
        $job_id="";
        unless (-f "$vcfname$sep_param$exe_condition.vcf")
        {
            if ($exe_condition eq "Default${sep_value}1")
        	{
        		$actual_exe_conditions="";
        		#print("DEBUG: A: Default exe conditions\n");
        	}
        	else
        	{
        		$actual_exe_conditions=join(" ",split("$sep_value",join(" ",split("$sep_param",$exe_condition))));
        		#print("DEBUG: A: Real exe_conditions\n");
        	}
            $job_id=submit_job_name($vcfname,"$qsub $variant_calling_sh $bamfiles $vcfname$sep_param$exe_condition.vcf $vcfname$sep_param${exe_condition}_platypus.log $actual_exe_conditions");
            print("\tSample $vcfname$sep_param$exe_condition.vcf variant calling submited with job_id $job_id\n");
        }
        else
        {
            print("\tSample $vcfname$sep_param$exe_condition.vcf variant calling already present. Skipping\n");
        }
        for my $filtering_condition (@filtering_conditions)
        {
            my $condition="$exe_condition$sep_param$filtering_condition";
            unless (-f "$vcfname$sep_param$condition.vcf")
            {
                my $filtering_command="$vcf_filt_exe ";
                $filtering_command.=join(" ",split("$sep_value",join(" ",split("$sep_param",$filtering_condition))));
                if($job_id eq "")
                {
                    #No dependency
                    $fjob_id=submit_job_name($vcfname,"$qsub_noparallel --job-name=$vcfname$sep_param${condition}_filtering $filtering_command -i $vcfname$sep_param$exe_condition.vcf -o $vcfname$sep_param$condition.vcf");
                    print("\tFiltering $vcfname$sep_param$exe_condition to generate $vcfname$sep_param$condition.vcf in job $fjob_id\n"); 
                }
                else
                {
                    $fjob_id=submit_job_name($vcfname,"$qsub_noparallel --dependency=afterok:$job_id --job-name=$vcfname$sep_param${condition}_filtering $filtering_command -i $vcfname$sep_param$exe_condition.vcf -o $vcfname$sep_param$condition.vcf");
                    print("\tFiltering $vcfname$sep_param$exe_condition to generate $vcfname$sep_param$condition.vcf in job $fjob_id\n"); 
                }
            }
            else
            {
                print("\tSample $vcfname$sep_param$condition.vcf already present. Skipping\n");
            }
        }
    }    
    for (my $j=$i+1; $j<scalar(@afiles);++$j)
    {
        $job_id="";
        $bamfiles="$normalfile,".$afiles[$i].",".$afiles[$j];
        $vcfname=basename($normalfile)."_".basename($afiles[$i])."_".basename($afiles[$j]);
        $vcfname=~s/.bam//g;
        unless (-f "$vcfname.vcf") ##Variant calling for each NAB file
        {
            $job_id=submit_job_name($vcfname,"$qsub $variant_calling_sh $bamfiles $vcfname.vcf ${vcfname}_platypus.log");
            print("\tSample $vcfname variant calling submited with job_id $job_id\n");
        }
        else
        {
            print("\tSample $vcfname variant calling already present. Skipping\n");
        }
    }
}

my $deps="";

foreach my $name (keys %cases)
{

    ($normal, $a, $b)=@{$cases{$name}};
    $normal=basename($normal);
    $normal=~s/.bam//;
    $a=basename($a);
    $a=~s/.bam//;
    $b=basename($b);
    $b=~s/.bam//;
    my @tempdeps;
    if(exists $job_ids{$normal})
    {
        push(@tempdeps,join(":",@{$job_ids{$normal}}));
    }
    if(exists $job_ids{$a})
    {
        push(@tempdeps,join(":",@{$job_ids{$a}}));
    }
    if(exists $job_ids{$b})
    {
        push(@tempdeps,join(":",@{$job_ids{$b}}));
    }
    if(exists $job_ids{"${normal}_${a}_${b}"}) ##
    {
        push(@tempdeps,join(":",@{$job_ids{"${normal}_${a}_${b}"}}));
    }
    if(exists $job_ids{"${normal}_${b}_${a}"})
    {
        push(@tempdeps,join(":",@{$job_ids{"${normal}_${b}_${a}"}}));
    }

    if (scalar @tempdeps > 0)
    {
        $deps="--dependency=afterok:".join(":",@tempdeps);
    }
    else
    {
        $deps="";
    }
    if(-f $onfile2)
    {   
        $job_id=submit_job_name("comp","$qsub $deps $helper_sh $helper_pl --output_folder $output_dir/$name -n $normal -a $a -b $b -e $oefile -f $offile --NABfilt_cond_inputfile $onfile --NABfilt_cond_inputfile2 $onfile2 --covaltB_cond_inputfile $ocfile -o $name.csv --n_cores $n_cores");
    }
    else
    {
        $job_id=submit_job_name("comp","$qsub $deps $helper_sh $helper_pl --output_folder $output_dir/$name -n $normal -a $a -b $b -e $oefile -f $offile --NABfilt_cond_inputfile $onfile --covaltB_cond_inputfile $ocfile -o $name.csv --n_cores $n_cores");
    }
    print("The filtering and analysis of the vcf files is being conducted with the job_id $job_id\n");
    
    $job_id=submit_job_name("tstv","$qsub_noparallel --dependency=afterok:$job_id $tstv_sh $output_dir/$name");
    print("Ts/Tv statistics for $name are being calculated in the job_id $job_id\n");
   
}

###Annovar and Ts/Tv for the main directory
#Depends on variant calling
my @temp;
foreach my $key (keys %job_ids)
{
    unless ($key =~ /tstv/) #All variant callings
    {
        push(@temp,@{$job_ids{$key}});
    }
}
if (scalar @temp > 0)
{
    $deps="--dependency=afterok:".join(":",@temp);
}
else
{
    $deps="";
}


#$job_id=submit_job_name("tstv","$qsub $deps $helper_sh $n_cores"); ##Annovar ##I think this is wrong

#$job_id=submit_job_name("tstv","$qsub_noparallel --dependency=afterok:$job_id $tstv_sh $output_dir"); #TsTv ##I think this is wrong

#print("Common Ts/Tv statistics are being calculated in the job_id $job_id\n"); ##What are these?

$deps=join(":",@{$job_ids{"tstv"}});

###EDIT: I definetly need to modify the postanalyzer, so that I can get all the data back.

### I will need to bring back the postanalyzer step. However, I do need to check other things first, this starting to be too complicated

#$job_id=submit_job_name("","$qsub_noparallel --dependency=afterok:$deps $postanalyzer_exe $output_dir $ifile $oefile $offile $onfile $onfile2 $ocfile");
#print("PostHeterAnalyzer job submitted with job_id $job_id\n");
#print("Jobs submitted!\n");

exit;

###################################################################################
###FUNCTIONS
###################################################################################
sub submit_job
{
    my ($command)=@_;
    my $job_id="";
    $job_id=`$command | $sed`;
    chomp($job_id);
    while ($job_id eq "")
    {
	    print("WARNING: Job submission has failed, trying again\n");
        sleep($sleep);
        $job_id=`$command | $sed`;
   		chomp($job_id);
    }
    #$job_ids{$job_id}=1;
    return $job_id;
}

sub submit_job_name
{
    my ($name,$command)=@_;
    my $job_id="";
    $job_id=`$command | $sed`;
    chomp($job_id);
    while ($job_id eq "")
    {
	    print("WARNING: Job submission has failed, trying again\n");
        sleep($sleep);
        $job_id=`$command | $sed`;
   		chomp($job_id);
    }
    if (exists $job_ids{$name})
    {
        push(@{$job_ids{$name}},$job_id);
    }
    else
    {
        $job_ids{$name}=[$job_id];
    }

    #print("DEBUG: job name $name, id $job_id\n");

    return $job_id;
}
sub compile_dependencies
{
    my $option="";
    if (scalar(@_) ==1)
    {
        ($option)=@_;
    }
    my $dependencies="";
    my $keep=0;
    if ($option =~/keep/i)
    {
        $keep=1;    
    }
    for my $key (keys %job_ids)
    {
        $dependencies="$dependencies$key:";
        $keep == 0 and delete($job_ids{$key});
    }
    chop($dependencies);
    return $dependencies;
}

sub wait_for_jobs
{
    my $awk='awk \'BEGIN{FS=" ";var=1}{if ($5 != "C"){var=0}}END{print var}\'';
    while (scalar keys %job_ids != 0)##Check 
    {
    	sleep($sleep); 
    	print("\tPending jobs ",join(",",keys %job_ids),"\n");
    	foreach my $id (keys %job_ids)
	    {
		    my $status=`$qstat $id | tail -n 1 | $awk`;
        	#print("DEBUG: Status job id $id : $status\n");
        	if($status eq 1) ## 1 means that the job has finished
        	{
            	delete($job_ids{$id});
        	}
	    }
    }
}

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
        @{ $ref_statistics }=(0,$n_common); ##Stats= proportion of selected reads in the reference, number of selected variants

    }
    else
    {
	    @{ $ref_statistics }=($n_common/$n2,$n_common); ##Stats= proportion of selected reads in the reference, number of selected variants
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
         @{ $ref_statistics }=(0,$n_selvariants); 
    }
    else
    {   	
	    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
	}
    return (\%variants1,\%variants2);
}

#Combine the variants of two samples. Calculate the union of the two samples
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
    @{ $ref_statistics }=($n_selvariants/($n_selvariants+$n_filtvariants),$n_selvariants); ##Stats= proportion of selected reads in the reference, number of selected variants
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
	for (my $i=1; $i<scalar @{$array};$i++) ###The first line is the header
	{
		$key=${$array}[$i];
		$key=~s/^(.*?)\t(.*?)\t.*/$1$OFS$2/;
		$hash{$key}=1;
		#print("DEBUG: New variant being hashed $key\n");
	}
	return \%hash;
}
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
#($ref_common_variantsAfilt,"Afilt$sep_param${condition}.vcf","A.vcf","##Comment for the header");
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

