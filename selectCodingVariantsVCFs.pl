#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;

my $usage="Usage: $0 [options] -i input_dir -o output_dir -f variantsfile\nThe variants file is the 1col result of tabulate_annovar_angelo.pl, filtered for the desired variants\n";
######################################################

######################################################
##MAIN
######################################################
##Configuration variables
my $IFS="\t";


##Input variables
my $inputDir="";
my $outputDir="";
my $variantsFile="";
my $help;

##Getopt
######################################################
(! GetOptions(
    'variants_file|f=s' => \$variantsFile,
    'output_dir|o=s' => \$outputDir,
    'input_dir|i=s' => \$inputDir,
    'help|h' => \$help,
                )) or (($outputDir eq "") || ($inputDir eq "")  || (! -s $variantsFile) || $help) and die $usage;


#Parse variantsFile
open(my $IFILE, $variantsFile);
my @variantsContent=<$IFILE>;
close($IFILE);

my @row;
my $case;
my $shared;
my $varid;
my $gene;
my $kind;
my $type;
my %variantdata;

for (my $i=1; $i<scalar @variantsContent; ++$i){ ##Skips the header
    @row=split("$IFS",$variantsContent[$i]);
    scalar @row != 8 and die "The line $i of the variant file, $variantsFile does not conform with the expected format\n";
    $case=$row[0];
    $shared=$row[1];
    $varid="$row[2]_$row[3]";
    $gene=$row[5];
    $kind=$row[6];
    $type=$row[7];

    if (exists $variantdata{$case}){
        if(exists $variantdata{$case}{$shared}){
            if(exists $variantdata{$case}{$shared}{$varid}){
                die "Error, several entries for the same variant";
            }else{
                $variantdata{$case}{$shared}{$varid}="$gene, $kind, $type";
            }
        }
        else{
            $variantdata{$case}{$shared}={$varid => "$gene, $kind, $type"};
        }
    }else{
        $variantdata{$case}={$shared => {$varid => "$gene, $kind, $type"}};
    }

}

##Main loop, directories
opendir(my $IDIR, $inputDir) or die "can't opendir $inputDir: $!";
my @dirs = grep { /^[^.]/ && -d "$inputDir/$_" } readdir($IDIR);
closedir $IDIR;

my $FILE;

foreach my $dir (@dirs){
    
    mkdir "$outputDir/$dir";
   
    #PrivateA
    filterVCF($dir,"A_private.vcf","A");

    #PrivateB
    filterVCF($dir,"B_private.vcf","B");
    
    #CommonA
    filterVCF($dir,"A_common.vcf","Common");
    
    #CommonB
    filterVCF($dir,"B_common.vcf","Common");

}    

sub filterVCF{
    my ($dir,$filename,$shared)=@_;
    open(my $FILE,"$inputDir/$dir/$filename");
    open(my $OFILE, ">$outputDir/$dir/$filename");
    my @ifile=<$FILE>;
    close($FILE);

    my @tempcontent;
    my $varid;

    foreach my $line (@ifile){
        if($line=~/^#/){
            print($OFILE $line);
        }else{
            @tempcontent=split($IFS,$line);
            $varid="$tempcontent[0]_$tempcontent[1]";
            if(exists $variantdata{$dir}{$shared}{$varid}){
                print($OFILE $line);
                #print("\nDEBUG: included line $line, varid: $varid");
            }else{
                #print("\nDEBUG: rejected line $line, varid: $varid");
            }
        }
    }
    close($OFILE);
    
}
