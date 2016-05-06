#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;
#use Parallel::Loops;

##Configuration variables
######################################################

our $sep_param="#";
our $sep_value=":";
our $OFS=",";
our $FS=",";
our $variant_caller="platypus";
our $variant_calling_sh;
our $helper_sh="vcfFilteringTableV2_analyser_helper.sh";
our $helper_pl="vcfFilteringTableV2_analyser.pl";
our $n_cores=1;
our $qsub="qsub -e e_logs/ -o o_logs/ -q shortq -l nodes=1:ppn=";

######################################################

##IO Variables
######################################################
my $execond_inputfile="";
my $filtercond_inputfile="";
my $NABfiltercond_inputfile="",
my $output_dir="vcf_outputdir";
my $output_file="";
my $normal_bam="";
my $sample1_bam="";
my $sample2_bam="";

#Flags
my $help;
my $usage="Usage: $0 [options] -o output_file --normal_bamfile bamfile_normal_sample --sample_A_bamfile bamfile_A_sample --sample_B_bamfile bamfile_B_sample\n\n\nOptions:\n--------\n\t-e/--exec_cond_inputfile : input file for execution parameters and options\n\t-f/--filt_cond_inputfile : input file for execution parameters and options\n\t--NABfilt_cond_inputfile : input file for the filtering options of the NAB sample\n\t--output_dir : output directory for vcf files\n\t--n_cores : number of cores to execute some steps in parallel\n\t\n\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
######################################################
(! GetOptions(
        'exec_cond_inputfile|e=s' => \$execond_inputfile,
	'filt_cond_inputfile|f=s' => \$filtercond_inputfile,
	'NABfilt_cond_inputfile=s' => \$NABfiltercond_inputfile,
        'output_dir=s' => \$output_dir,
	'output_file|o=s' => \$output_file,
	'normal_bamfile=s' => \$normal_bam,
	'sample_A_bamfile=s' => \$sample1_bam,
	'sample_B_bamfile=s' => \$sample2_bam,
        'n_cores=i' => \$n_cores,
        'help|h' => \$help,
                )) or (($output_file eq "") || ($normal_bam eq "")  || ($sample1_bam eq "") || ($sample2_bam eq "")  || $help) and die $usage;

$qsub.=$n_cores;

##Input file parsing and directory creation
######################################################

my @exe_parameters=("input");
my @exe_param_values=([("")]);

##Input files

if ($execond_inputfile ne "")
{
	@exe_parameters=();
	parse_parameters_values($execond_inputfile,\@exe_parameters,\@exe_param_values);
}

##BAMs
if ( ! -f $normal_bam)
{
	die "The BAM file $normal_bam is not accesible. Please, check your input options.\n";
}
else
{
	$normal_bam=Cwd::abs_path($normal_bam);
}

if ( ! -f $sample1_bam)
{
	die "The BAM file $sample1_bam is not accesible. Please, check your input options.\n";
}
else
{
	$sample1_bam=Cwd::abs_path($sample1_bam);
}

if ( ! -f $sample2_bam)
{
	die "The BAM file $sample2_bam is not accesible. Please, check your input options.\n";
}
else
{
	$sample2_bam=Cwd::abs_path($sample2_bam);
}

mkdir $output_dir;
my $original_dir=dirname(Cwd::abs_path($0));
my $oefile=Cwd::abs_path($execond_inputfile);
my $offile=Cwd::abs_path($filtercond_inputfile);
my $onfile=Cwd::abs_path($NABfiltercond_inputfile);
chdir $output_dir or die "The output directory $output_dir is not accesible";

#my $vcf_filt_exe="vcf_filtering.pl";
#
#if(`which vcf_filtering.pl 2>/dev/null` eq "")
#{
#	if (-f "$original_dir/$vcf_filt_exe")
#	{
#		$vcf_filt_exe="$original_dir/$vcf_filt_exe";
#	}
#	else 
#	{
#		die "The executable vcf_filtering.pl is neither in your PATH nor in the directory $original_dir. This script needs to locate it to continue\n";
#	}
#}

if (-f "$original_dir/$variant_caller.sh")
{
	$variant_calling_sh="$original_dir/$variant_caller.sh";
}
else
{
	die "Error, the sh file for the variant caller $variant_caller is not located in $original_dir, named $original_dir/$variant_caller.sh. Please, fix this in order to use this script\n"
}

if (-f "$original_dir/$helper_sh" && -f "$original_dir/$helper_pl")
{
    $helper_sh="$original_dir/$helper_sh";
    $helper_pl="$original_dir/$helper_pl";
}
else
{
	die "Error, the sh file for the secondary analysis of the data, $helper_sh is not located in $original_dir, named $original_dir/$helper_pl. Please, fid this in order to use this script\n"
}

## Main conditions loop
#######################

my $name_condition;
my @exe_conditions;
my $bamfiles;

##Generate all combinations of proposed values for execution
combs(0,"",\@exe_parameters,\@exe_param_values,\@exe_conditions);

my %job_ids;
my $exe_condition;

mkdir("e_logs");
mkdir("o_logs");

### Default calling parameters. This may be changed in the future
#########################################################

print("Unfiltered variant calling detection/execution:\n");	
if(! -f "N.vcf")
{
	$bamfiles=$normal_bam;
	$exe_condition="";
	my $job_id=`$qsub $variant_calling_sh -F "$bamfiles $exe_condition N.vcf N_platypus.log" | sed "s/.master.cm.cluster//"`;
   		chomp($job_id);
    	$job_ids{$job_id}=1;
	print("\tNormal tissue variant calling submited with job_id $job_id\n");
}
else
{
	print("\tNormal tissue variant calling already present, skipping it\n");
}

if(! -f "A.vcf")
{
	$bamfiles=$sample1_bam;
	$exe_condition="";
	my $job_id=`$qsub $variant_calling_sh -F "$bamfiles $exe_condition A.vcf A_platypus.log" | sed "s/.master.cm.cluster//"`;
    	chomp($job_id);
    	$job_ids{$job_id}=1;
	print("\tSample A variant calling submited with job_id $job_id\n");

}
else
{
	print("\tSample A variant calling already present, skipping it\n");
}

if(! -f "B.vcf")
{
	$bamfiles=$sample2_bam;
	$exe_condition="";
	my $job_id=`$qsub $variant_calling_sh -F "$bamfiles $exe_condition B.vcf B_platypus.log" | sed "s/.master.cm.cluster//"`;
    	chomp($job_id);
    	$job_ids{$job_id}=1;
	print("\tSample B variant calling submited with job_id $job_id\n");

}
else
{
	print("\tSample B variant calling already present, skipping it\n");
}

if(! -f "NAB.vcf")
{
	$bamfiles="$normal_bam,$sample1_bam,$sample2_bam";
	$exe_condition="";
	my $job_id=`$qsub $variant_calling_sh -F "$bamfiles $exe_condition NAB.vcf NAB_platypus.log" | sed "s/.master.cm.cluster//"`;
    	chomp($job_id);
    	$job_ids{$job_id}=1;
	print("\tSample NAB variant calling submited with job_id $job_id\n");

}
else
{
	print("\tNAB variant calling already present, skipping it\n");
}

### Calling with filters (Only cancer samples so far)
##############################################################
my $job_id;
	my $actual_exe_conditions;
print("Filtered variant calling detection/execution:\n");	
foreach my $exe_condition (@exe_conditions)
{
	if (! -f "A$sep_param$exe_condition.vcf")
	{
		$bamfiles=$sample1_bam;
        		if ($exe_condition eq "Default${sep_value}1")
        		{
            		$actual_exe_conditions="";
            		#print("DEBUG: A: Default exe conditions\n");
       	 	}
        		else
        		{
            		$actual_exe_conditions=$exe_condition;
            		#print("DEBUG: A: Real exe_conditions\n");
        		}
        		#print("DEBUG: qsub -e e_logs/ -o o_logs/ -q shortq $variant_calling_sh -F \"$bamfiles $actual_exe_conditions A$sep_param$exe_condition.vcf A$sep_param${exe_conditions}_platypus.log\" | sed \"s/.master.cm.cluster//\"");
		$job_id=`$qsub $variant_calling_sh -F "$bamfiles $actual_exe_conditions A$sep_param$exe_condition.vcf A$sep_param${exe_condition}_platypus.log" | sed "s/.master.cm.cluster//"`;
		chomp($job_id);
        		$job_ids{$job_id}=1;
		print("\tSample A variant calling for conditions $exe_condition submited with job_id $job_id\n");

	}
	else
	{
		print("\tSample A variant calling for conditions $exe_condition already present, skipping it\n");
	}

	if (! -f "B$sep_param$exe_condition.vcf")
	{
		$bamfiles=$sample2_bam;
        		if ($exe_condition eq "Default${sep_value}1")
        		{
            		$actual_exe_conditions="";
            		##print("DEBUG: B: Default exe conditions\n");
        		}
        		else
        		{
           	 		$actual_exe_conditions=$exe_condition;
            		#print("DEBUG: B: Real exe_conditions\n");
        		}
		$job_id=`$qsub $variant_calling_sh -F "$bamfiles $actual_exe_conditions B$sep_param$exe_condition.vcf B$sep_param${exe_condition}_platypus.log" | sed "s/.master.cm.cluster//"`;
    		chomp($job_id);	
        		$job_ids{$job_id}=1;
		print("\tSample B variant calling for conditions $exe_condition submited with job_id $job_id\n");

	}
	else
	{
		print("\tSample B variant calling for conditions $exe_condition already present, skipping it\n");

	}

} 

while (scalar keys %job_ids != 0)##Check 
{
    	sleep(60); 
    	print("\tPending jobs ",join(",",keys %job_ids),"\n");
	foreach my $id (keys %job_ids)
	{
		my $status=system "qstat $id >/dev/null 2>&1";
        		#print("DEBUG: Status job id $id : $status\n");
        		if($status!=0)
        		{
            		delete($job_ids{$id});
        		}
	}
	
}

$job_id=`$qsub $helper_sh -F "$helper_pl -e $oefile -f $offile --NABfilt_cond_inputfile $onfile -o $output_file --original_directory $original_dir" | sed "s/.master.cm.cluster//"`;
chomp($job_id);

print "The filtering and analysis of the vcf files is being conducted with the job_id $job_id. Now this script will exit\n";
exit;

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

