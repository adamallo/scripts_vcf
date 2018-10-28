#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Bio::DB::HTS::Tabix;

#use Cwd;
#use File::Basename;
#use Env;
#use Parallel::Loops;

##Configuration variables
######################################################
my $IFS="\t";
my $OFS="\t";
my $AFS=","; ##Allele separation
######################################################

##IO Variables
######################################################
my $pop_file="";
my $input_file="";
my $output_file="";

#Flags
my $help;
my $usage="Usage: $0 -p populationDataVCFWith_Tabix -i input.vcf -o output.tsv\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
######################################################
(! GetOptions(
	'population_file|p=s' => \$pop_file,
    'input_file|i=s' => \$input_file,
    'output_file|o=s' => \$output_file,
    'help|h' => \$help,
                )) or (($output_file eq "") || ($input_file eq "")  || ($pop_file eq "") || $help) and die $usage;


##Input file parsing and directory creation
######################################################

-s $pop_file && -s "${pop_file}.tbi" or -s die "The file with population data, $pop_file cannot be read or does not have the proper tabix format";

my $tabix = Bio::DB::HTS::Tabix->new( filename =>$pop_file );

open(my $FILEHANDLE, $input_file) or die "The input file $input_file does not exist or cannot be read";
my @input_data=<$FILEHANDLE>;
close($FILEHANDLE);

my $i=0;
my $chr="";
my $nstart=0;
my $tstart=0;
my $nend="";
my $ref="";
my $tref="";
my $ialt="";
my $talt="";
my $tfilt="";
my $taf=0;
my $tabix_iter;
my %pop_data;

open(my $OUTPUT_FILE, ">$output_file") or die "The output file $output_file cannot be opened";
print($OUTPUT_FILE "#CHROM\tPOS\tREF\tALT\tFILTER\tAF\n");

for ($i=0; $i<scalar @input_data; ++$i) 
{
    if($input_data[$i]=~m/^[^#]/)
    {
        my @linecontent=split($IFS,$input_data[$i]);
        $chr=$linecontent[0];
        $nstart=$linecontent[1];
        $ref=$linecontent[3];
        my @alt=split($AFS,$linecontent[4]);
        
        #print("DEBUG: $chr $nstart $ref @alt\n");
        #print("DEBUG: "."$chr:$nstart-".($nstart+1));
        $tabix_iter=$tabix->query("$chr:$nstart-".($nstart+1));
        
        if(defined $tabix_iter)
        { 
            %pop_data=();
            while(my $line=$tabix_iter->next)
            {
                #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  S1 ... SN
                $line =~ s/^[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t[^\t]+(AF=[^;]+).*$/$1\t$3\t$2\t$4\t$5/;
                #print("DEBUG: $line\n");
                ($tstart,$talt,$tref,$tfilt,$taf)=split($IFS,$line);
                #print("DEBUG: $tstart, $talt, $tref, $tfilt, $taf\n");
                $pop_data{"${tstart}_$talt"}=[$tref,$tfilt,$taf];
            }
           
        }
        else
        {
            warn "There is no data for $chr:$nstart-".($nstart+1);
        }
        foreach $ialt (@alt) 
        {
            if(exists $pop_data{"${nstart}_$ialt"}){
               print($OUTPUT_FILE "$chr$OFS$nstart$OFS$ref$OFS$ialt$OFS".$pop_data{"${nstart}_$ialt"}[1].$OFS.$pop_data{"${nstart}_$ialt"}[2]."\n"); 
            }
            else
            {
               print($OUTPUT_FILE "$chr$OFS$nstart$OFS$ref$OFS$ialt${OFS}NA${OFS}NA\n"); 
            }
        }
    }
    
}

close($OUTPUT_FILE);
