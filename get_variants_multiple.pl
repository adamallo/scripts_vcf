use warnings;
use strict;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use File::Path;

#IO
my $infile="";
my $indir="";
my $outdir="";
my $covB="";
my $PAF="";

#Configuration variables
my $OFS="\t";
my $IFS="\t";

#Flags
my $help;
my $usage="Usage: $0 -f file_list -d input_dir -o output_dir [options]\n-f file list with same format as for HeterAnalyzer_loop_multiple.sh\n-d output dir indicated to HeterAnalyzer_loop_multiple.sh\n-o directory in which this script will generate a vcf file per sample and an additional statistics file\n\nOptions:\n--------\n\t--covB: suffix to select the covB results\n\t--PAF: suffix to select the PAF results\n\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
######################################################
(! GetOptions(
    'files|f=s' => \$infile,
    'input_dir|d=s' => \$indir,
    'output_dir|o=s' => \$outdir,
    'covB|b=s' =>\$covB,
    'PAF|p=s' =>\$PAF,
    'help|h' => \$help,
                )) or (($outdir eq "") || ($indir eq "") || ($infile eq "") || $help) and die $usage;


##Main data structures
my @dicts; #Array of dict arrays, one per sample, with information on all their variants
my %vars; #Dict of variants, with information on the number of samples that carry them
my %samples; #Dict of samples (SXX format), with refs to arrays with the number of clonal, subclonal, and private mutations
my $nsamples;
my @samplename; #Array of sample names
my @ovcfs; #Array of the original vcf filenames for each sample
my @folders;
my $ndigits;

##Parsing infile
mopen(my $FH_INFILE, $infile);
my @content=<$FH_INFILE>;
close($FH_INFILE);
my ($name,$normal)=splice(@content,0,2);
chomp($name);
chomp($normal);
$nsamples=scalar @content;
$ndigits=length $nsamples;

print("Parsing input files from all comparisons of $nsamples samples in the folder $indir\n");

#Parsing vcf files

#Initializing dictionaries
for (my $i=0; $i<$nsamples; ++$i)
{
    $dicts[$i]={};
}

##Main parsing loop
for (my $i=0; $i<$nsamples-1; ++$i)
{
    my $Asample=sprintf("S%0*d",$ndigits,$i);
    my $ovcf;
    
    for (my $j=$i+1; $j<$nsamples; ++$j)
    {
        #Identify A and B samples
        my $Bsample=sprintf("S%0*d",$ndigits,$j);
        my $folder=join("_",$name,$Asample,$Bsample);
        push(@folders,$folder);
        
        #Look for the files that containt the A and B variants we want
        open(my $vcfdictFH, "$indir/$folder/vcfdict.csv") or die "ERROR: $indir/$folder/vcfdict.csv cannot be opened\n";
        my @vcfdictcontent=<$vcfdictFH>;
        close($vcfdictFH);

        print("\tFolder $folder\n");
    
        my @afiles=get_files(\@vcfdictcontent,"Afilt${covB}NAB${PAF}#.*different.vcf");
        my @bfiles=get_files(\@vcfdictcontent,"Bfilt${covB}NAB${PAF}#.*different.vcf");
        my @commonfiles=get_files(\@vcfdictcontent,"^filt${covB}NABU${PAF}#.*common.vcf");

        scalar  @afiles != 1 || scalar @bfiles != 1 || scalar @commonfiles != 1 and die "Error detecting input vcf files in folder $folder\n";
    
        add_vcfvariants($dicts[$i],\@afiles);
        add_vcfvariants($dicts[$j],\@bfiles);
       
        #I don't have a subroutine for two but I could 
        my @commonvariants=parse_vcf($commonfiles[0]);
        foreach my $var (@commonvariants)
        {
            $dicts[$i]->{$var}+=1;
            $dicts[$j]->{$var}+=1;
        }

    
        #Getting the name of Asample and the original vcf file 
        if($j==$i+1) ##In the inside loop because we need to parse one of the specific vcfdict files 
        {
            ($ovcfs[$i],$samplename[$i])=get_vcf_samplename(",.*${Asample}#",\@vcfdictcontent);##TODO:I know this one is originally as an absolute path. This is a mess and should be solved

            print("\t\t$Asample name $samplename[$i], original vcf file $ovcfs[$i]\n");
            if($j==$nsamples-1)
            {
               ($ovcfs[$j],$samplename[$j])= get_vcf_samplename(",.*${Bsample}#",\@vcfdictcontent);##TODO:I know this one is originally as an absolute path. This is a mess and should be solved

            print("\t\t$Bsample name $samplename[$j], original vcf file $ovcfs[$j]\n");
            }
        }
    }
}

print("Done\n");

print("Generating statistics and writting output files...");

##Generate statistics
##Dictionary of variants, with number of samples that bear them
for (my $isample=0; $isample<$nsamples; ++$isample)
{
    foreach my $var (keys %{$dicts[$isample]})
    {
        $vars{$var}+=1;
    }
    
} #foreach my $sample

##Dictionary of samples, with the number of private, shared, clonal mutation
for (my $isample=0; $isample<$nsamples; ++$isample)
{
    my $Asample=sprintf("S%0*d",$ndigits,$isample);
    $samples{$Asample}=[0,0,0];
    foreach my $var (keys %{$dicts[$isample]})
    {
        if($vars{$var}>1)
        {
            if($vars{$var}==$nsamples)
            {
                $samples{$Asample}->[0]+=1;
            }
            else
            {
                $samples{$Asample}->[1]+=1;
            }
        }
        else
        {
            $samples{$Asample}->[2]+=1;
        }
    }
    
} #foreach my $sample

##Outputs
mkpath("$outdir/vcfs");
mkpath("$outdir/stats");
# 1 samples file with the number of clonal, subclonal, and private variants
mopen(my $OUT_SAMPLE_STATS, ">$outdir/stats/sample.csv");
print($OUT_SAMPLE_STATS "Sample,Clonal,Subclonal,Private\n");
for (my $isample=0; $isample<$nsamples; ++$isample)
{
    my $sample=sprintf("S%0*d",$ndigits,$isample);
    my $outfile="$outdir/vcfs/${sample}_final.vcf";

    # 1 vcf file per sample
    write_variant_vcf($dicts[$isample],$outfile,$ovcfs[$isample],"#Variants filtered using HeterAnalyzer_multiple.pl");

    # 1 line per sample summary
    print($OUT_SAMPLE_STATS join(",", $samplename[$isample], @{$samples{$sample}})."\n");
}
close($OUT_SAMPLE_STATS);

#TODO: # 1 csv file per sample with the important variant information, with number of reads of reference and alt and same for normal?

# 1 variants file with the number of samples that contain it and a classification of clonal, subclonal, and private
mopen(my $OUT_VAR_STATS, ">$outdir/stats/var.csv");
print($OUT_VAR_STATS "CHROM,POS,ID,REF,ALT,NSAMPLES,KIND\n");
foreach my $var (keys %vars)
{
    my @temp=split($OFS,$var);
    my $kind="Private";

    if($vars{$var}>1)
    {
        $kind=$vars{$var}==$nsamples?"Clonal":"Subclonal";
    }
    print($OUT_VAR_STATS join(",",@temp[0..1],$var,@temp[2..3],$vars{$var},$kind)."\n");
}
close($OUT_VAR_STATS);

print(" Done\n");

exit;

##FUNCTIONS
sub get_files
{
    my ($refcontent,$regex,$dir)=@_;
    my @output;
    my @temp;

    if(defined $dir)
    {
        $dir=$dir."/";
    }
    else
    {
        $dir="";
    }

    for my $line (@{$refcontent})
    {
        chomp($line);
        #print("DEBUG: $line\n");
        if ($line =~ m/$regex/)
        {
            @temp=split(",",$line);
            push(@output,$dir.$temp[1]);
           #print("\tDEBUG: pushing $temp[1] in common files\n");
        }
    }
    
    return @output;
}

# Parses the variants in vcf files pointed by the array @{$filesref} and adds them to the hash %{$hashref}
# It depends on parse_vcf
# The value indicates the number of comparisons in which that variant is present. I could use this as a
# measure of support somehow.
# ########################################################################################################

sub add_vcfvariants
{
    my ($hashref,$filesref)=@_;

    foreach my $file (@{$filesref})
    {
        my @variants=parse_vcf($file);
        foreach my $var (@variants)
        {
            $hashref->{$var}+=1;
        }
    }

}

# This function uses the parsed info from a vcfdict file to obtain the filename of
# the original vcf file for a sample and the name of that sample, returning them
sub get_vcf_samplename
{
    my ($regex,$refvcfdictcontent,$dir)=@_;
    my @retdata;
    my @temp=get_files($refvcfdictcontent,$regex,$dir);
    scalar @temp != 1 and die "ERROR: detected more than one candidate vcf file for sample $regex\n";
    push(@retdata,$temp[0]);
    @temp=get_sample_names_vcf($retdata[0]);
    scalar @temp != 1 and die "ERROR: detected more than one sample in the vcf file $retdata[0]\n";
    $temp[0]=~s/H_SL-DCIS-//g; ##Simplification to make them easier to read #TODO: HARDCODED
    push(@retdata,$temp[0]);
    return @retdata;
}


##THESE SUBRUTINES SHOULD PROBABLY BE PART OF MY OWN MODULE/LIBRARY

# Parse vcf generating variant keys in the form CHROM$OFSPOS$OFSREF$OFSALT
##########################################################################
sub parse_vcf
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file cannot be opened";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my @keys;
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
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            chomp($key);
            push(@keys,$key);
        }
    }
    return @keys;
}

# Writes a vcf with the variables contained in a hash selected from another VCF file
# ##################################################################################

sub write_variant_vcf
{
    my ($ref_hash,$filename,$vcf,$comment)=@_;
    open(my $OFILE, ">$filename") or die "ERROR: The file $filename could not be opened for writting";
    open(my $IFILE, "$vcf");
    my @icontent= <$IFILE>;
    close($IFILE);
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
            chomp($key);
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
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

# Parses the sample names from a VCF and returns an array with them

sub get_sample_names_vcf
{
    my ($vcfname)=@_;
    my @names;
    
    mopen(my $FT, $vcfname);
    my @content=<$FT>;
    close($FT);
    
    foreach my $line (@content)
    {
        if($line =~m/^#[^#]/)
        {
            chomp($line);
            @names=split($IFS,$line);
            splice(@names,0,9);
            last;
        }
    }
    return @names;
}

# Open with error message
###################################################################################
sub mopen
{
    open($_[0],$_[1]) or die "ERROR: impossible to open the file $_[1] in ".($_[1]=~m/^>/?"write":"read")."mode\n";
    return $_[0];
}

