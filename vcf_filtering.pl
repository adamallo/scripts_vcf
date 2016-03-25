#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

my $SNPSIFT_DIR=$ENV{'SNPSIFT_DIR'};
my ($qual,$AtoC,$AtoG,$AtoT,$CtoA,$CtoG,$CtoT,$GtoA,$GtoC,$GtoT,$TtoA,$TtoC,$TtoG,$min_cover); #Only some options to start with
my $input_file="";
my $output_file="";
my $help;
my $usage="Usage: $0 -i input_file -o output_file [options]\n\n\nOptions:\n--------\n\t-q/--qual : quality filter\n\t--atoc : filter out mutations from A to C\n\t--atog :  filter out mutations from A to G\n\t--atot : filter out mutations from A to T\n\t--ctoa :  filter out mutations from C to A\n\t--ctog : filter out mutations from C to G\n\t--ctot :  filter out mutations from C to T\n\t--gtoa : filter out mutations from G to A\n\t--gtoc :  filter out mutations from G to C\n\t--gtot : filter out mutations from G to T\n\t--ttoa :  filter out mutations from T to A\n\t--ttoc : filter out mutations from T to C\n\t--ttog : filter out mutations from T to G\n\t-m/--min_coverage :minimum coverage per locus\n\t\n\n";
my $exe;

########### SNPSIFT Util detection #########################
if (defined $SNPSIFT_DIR)
{
	if (-f $SNPSIFT_DIR."SnpSift.jar")
	{
		$exe=$SNPSIFT_DIR."SnpSift.jar";
	}
	elsif (-f "$SNPSIFT_DIR/SnpSift.jar")
	{
		$exe="$SNPSIFT_DIR/SnpSift.jar";
	}
}
elsif (-f "SnpSift.jar")
{
	$exe="SnpSift.jar";
}
else
{
	die "ERROR: SnpSift.jar not found. Please, specify its location using the environmental variable SNPSIFT_DIR\n";
}
############################################################

##################     ARGV parsing   ######################

((! GetOptions(
	'input|i=s' => \$input_file,
	'output|o=s' => \$output_file,
	'qual|q=i' => \$qual,
	'atoc' => \$AtoC,
	'atog' => \$AtoG,
	'atot' => \$AtoT,
	'ctoa' => \$CtoA,
	'ctog' => \$CtoG,
	'ctot' => \$CtoT,
	'gtoa' => \$GtoA,
	'gtoc' => \$GtoC,
	'gtot' => \$GtoT,
	'ttoa' => \$TtoA,
	'ttoc' => \$TtoC,
	'ttog' => \$TtoG,
	'min_cover|m=i' => \$min_cover, 
	'help|h' => \$help,
		)) or ((! -f $input_file) || ($output_file eq "") || $help)) and die $usage;

my ($qual_filter,$ctot_filter,$gtoa_filter,$min_cover_filter);

################### VCF scan ##############################

open(my $INPUTFILE,$input_file);
my $backup=$/;
$/="";
my $content=<$INPUTFILE>;
$/=$backup;
close($INPUTFILE);

my $platypus=0;
my $multisnv=0;

if ($content=~/platypus/i)
{
	$platypus=1;	
}
elsif ($content=~/multisnv/i)
{
	$multisnv=1;
}
else
{
	die "The input VCF file has not been recognized as generated by either platypus or multisnv. This script has been intended only for filtering VCFs generated by those pieces of software.\n";
}

###########################################################

################## Filter generation #####################

my $out_filter="";
my $variable="";

if ($qual)
{
	$out_filter.=" ( QUAL >= $qual ) &";
}
if ($AtoC)
{
	$out_filter.=" ( REF!=\'A\' || ( REF==\'A\' & ALT!=\'C\' ) ) &";
}
if ($AtoG)
{
	$out_filter.=" ( REF!=\'A\' || ( REF==\'A\' & ALT!=\'G\' ) ) &";
}
if ($AtoT)
{
	$out_filter.=" ( REF!=\'A\' || ( REF==\'A\' & ALT!=\'T\' ) ) &";
}
if ($CtoA)
{
	$out_filter.=" ( REF!=\'C\' || ( REF==\'C\' & ALT!=\'A\' ) ) &";
}
if ($CtoG)
{
	$out_filter.=" ( REF!=\'C\' || ( REF==\'C\' & ALT!=\'G\' ) ) &";
}
if ($CtoT)
{
	$out_filter.=" ( REF!=\'C\' || ( REF==\'C\' & ALT!=\'T\' ) ) &";
}
if ($GtoA)
{
	$out_filter.=" ( REF!=\'G\' || ( REF==\'G\' & ALT!=\'A\' ) ) &";
}
if ($GtoC)
{
	$out_filter.=" ( REF!=\'G\' || ( REF==\'G\' & ALT!=\'C\' ) ) &";
}
if ($GtoT)
{
	$out_filter.=" ( REF!=\'G\' || ( REF==\'G\' & ALT!=\'T\' ) ) &";
}
if ($TtoA)
{
	$out_filter.=" ( REF!=\'T\' || ( REF==\'T\' & ALT!=\'A\' ) ) &";
}
if ($TtoC)
{
	$out_filter.=" ( REF!=\'T\' || ( REF==\'T\' & ALT!=\'C\' ) ) &";
}
if ($TtoG)
{
	$out_filter.=" ( REF!=\'T\' || ( REF==\'T\' & ALT!=\'G\' ) ) &";
}
if ($min_cover)
{
	if ($platypus)
	{
		$variable="TC";
	}
	elsif ($multisnv)
	{
		$variable="DP"
	}
	$out_filter.=" ( $variable >= $min_cover) &";
}

$out_filter eq "" and die "ERROR: No filtering options have been specified\n\n$usage"; 

chop($out_filter); ##Removing extra " &"
chop($out_filter);

#########################################################

################      SNPsift execution      ############

substr($out_filter,0,1)=""; ##Removing extra space in the begining
system("cat $input_file | java -jar $exe filter \"$out_filter\" > $output_file");
exit;

#########################################################
