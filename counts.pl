###
#Perl script that filters out the content of one tsv file based on the presence in another one, using a reference column
###

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#Configuration variables
our $FS="\t";
our $OFS="\t";

#IO
my $usage="$0 --inputfile1 input_file1 -o ouput_file_prefix\n";
my $inputfile1="";
my $outputfile_prefix="";
my $help;

##Getopt
######################################################
(! GetOptions(
'inputfile1=s' => \$inputfile1,
'outputfile_prefix|o=s' => \$outputfile_prefix,
'help|h' => \$help,
)) or (($outputfile_prefix eq "") || $help ) and die $usage;

open(my $INPUT1,$inputfile1) or die "The input file $inputfile1 cannot be opened.\n $usage";

open(my $OUTPUT1,">${outputfile_prefix}_summary.tsv") or die "The output file ${outputfile_prefix}_summary.tsv cannot be opened.\n $usage";
open(my $OUTPUT2,">${outputfile_prefix}_bysample.tsv") or die "The output file ${outputfile_prefix}_bysample.tsv cannot be opened.\n $usage";

open(my $OUTPUT1GENE,">${outputfile_prefix}_summary_gene.tsv") or die "The output file ${outputfile_prefix}_summary_gene.tsv cannot be opened.\n $usage";
open(my $OUTPUT2GENE,">${outputfile_prefix}_bysample_gene.tsv") or die "The output file ${outputfile_prefix}_bysample_gene.tsv cannot be opened.\n $usage";

#Input
######
my @data1=<$INPUT1>;
close $INPUT1;

#Hashing
########

my %sumdata;
my %sample_dict;

my %sumdata_gene;
my %sample_dict_gene;

my $a;
my $b;
my $name;
my $gene;
my @temp;
my $kind;
my $tempn;
my $tempstring;

my $n=0;
my $samples="";
my $otherinfo="";

shift @data1; #Eliminating the header

foreach my $line (@data1)
{
    chomp($line);
    @temp=split($FS,$line);
    $a=$temp[0];
    $b=$temp[1];
    $kind=$temp[2];
    $name="$temp[3]_$temp[4]";
    $gene=$temp[6];
    $otherinfo="NA${OFS}NA";
    if(scalar @temp == 8)
    {
        $otherinfo="$temp[7]${OFS}NA";
    }
    elsif(scalar @temp == 9)
    {
        if($temp[8] eq "")
        {
            $temp[8]="NA";
        }
        $otherinfo="$temp[7]$OFS$temp[8]";
    }

    #print("DEBUG: $a, $b, $kind, $name, $gene, $otherinfo\n");
    if ($kind eq "A")
    {
        $n=1;
        $samples=$a;
    }
    elsif ($kind eq "B")
    {
        $n=1;
        $samples=$b;
    }
    elsif ($kind eq "common")
    {
        $n=2;
        $samples="$a$OFS$b";
    }
    else
    {
        die "Error, type $kind not compatible with this script\n";
    }
    unless (exists $sumdata{$name})
    {
        #print("DEBUG: New entry, $name, $n, $samples\n");
        $sumdata{$name}=[$n,$samples,"$gene$OFS$otherinfo"];
    }
    else 
    {
        $tempn=${$sumdata{$name}}[0];
        $tempstring=${$sumdata{$name}}[1];
        $tempn=$tempn+$n;
        $tempstring=join($OFS,$tempstring,$samples);
        #print("DEBUG: Updating entry, $name, ${$sumdata{$name}}[0] to $tempn, ${$sumdata{$name}}[1] to $tempstring\n");
        ${$sumdata{$name}}[0]=$tempn;
        ${$sumdata{$name}}[1]=$tempstring;
    }
    $sample_dict{$a}=1;
    $sample_dict{$b}=1;
    
    unless (exists $sumdata_gene{$gene})
    {
        #print("DEBUG: New entry, $gene, $n, $samples\n");
        $sumdata_gene{$gene}=[$n,$samples];
    }
    else 
    {
        $tempn=${$sumdata_gene{$gene}}[0];
        $tempstring=${$sumdata_gene{$gene}}[1];
        $tempn=$tempn+$n;
        $tempstring=join($OFS,$tempstring,$samples);
        #print("DEBUG: Updating entry, $gene, ${$sumdata_gene{$gene}}[0] to $tempn, ${$sumdata_gene{$gene}}[1] to $tempstring\n");
        ${$sumdata_gene{$gene}}[0]=$tempn;
        ${$sumdata_gene{$gene}}[1]=$tempstring;
    }
}


foreach my $key (keys %sumdata)
{
    print($OUTPUT1 join($OFS,$key,${$sumdata{$key}}[2],${$sumdata{$key}}[0]),"\n");
}

foreach my $key (keys %sumdata_gene)
{
    print($OUTPUT1GENE join($OFS,$key,${$sumdata_gene{$key}}[0]),"\n");
}

my @uniqsamples=keys %sample_dict;

#print("DEBUG: ",join($OFS,@uniqsamples),"\n");


my $tempout;
print($OUTPUT2 join($OFS,"Sample","Gene","Location","Type",@uniqsamples),"\n");
foreach my $key (keys %sumdata)
{
    $tempout="";
    foreach my $sample (@uniqsamples)
    {
        $n=0;
        foreach my $field (split($OFS,${$sumdata{$key}}[1]))
        {
            if ($field eq $sample)
            {
                $n+=1;
            }
        }
        $tempout.="$n$OFS";
    }
    chop($tempout);
    print($OUTPUT2 join($OFS,$key,${$sumdata{$key}}[2],$tempout),"\n");
}

print($OUTPUT2GENE join($OFS,"Sample",@uniqsamples),"\n");
foreach my $key (keys %sumdata_gene)
{
    $tempout="";
    foreach my $sample (@uniqsamples)
    {
        $n=0;
        foreach my $field (split($OFS,${$sumdata_gene{$key}}[1]))
        {
            if ($field eq $sample)
            {
                $n+=1;
            }
        }
        $tempout.="$n$OFS";
    }
    chop($tempout);
    print($OUTPUT2GENE join($OFS,$key,$tempout),"\n");
}

close $OUTPUT1;
close $OUTPUT2;
close $OUTPUT1GENE;
close $OUTPUT2GENE;
