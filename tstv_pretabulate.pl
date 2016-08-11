#!/bin/perl -w

use warnings;
use strict;

#Configuration parameters
#########################

our $sep_param="#";
our $sep_value=":";
our $OFS=",";
our $FS=",";

###CHECK ARGV
if(scalar @ARGV != 2)
{
    die "Usage: script inputfile outputfile\n";
}
my $inputfile=$ARGV[0];
my $outputfile=$ARGV[1];

open(my $IFILE, $inputfile) or die "Error reading the input file $inputfile\n";
my @content=<$IFILE>;
close($IFILE);

my %fullconds;
my %Aconds;
my %Bconds;
my %filtconds;
my %filtNconds;
my @temp;
my $cond;

foreach my $line (@content)
{
    chomp($line);
    @temp=split($FS,$line);
    if($temp[0]=~s/^filtNAB${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $cond=(split("${sep_param}NAB${sep_param}",$temp[0]))[0];
        if(defined $fullconds{$cond})
        {
            #print("DEBUG $cond, $fullconds{$cond}\n");
            ${$fullconds{$cond}}{$temp[0]}=$temp[1];
        }
        else
        {
            my %hash;
            $hash{$temp[0]}=$temp[1];
            $fullconds{$cond}=\%hash;
        }
    }
    elsif($temp[0]=~s/^A${sep_param}//)
    {
        $Aconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^B${sep_param}//)
    {
        $Bconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtN${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtNconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filt${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtconds{$temp[0]}=$temp[1];
    }
}

open(my $OFILE, ">$outputfile") or die "The file $outputfile can't be opened, check your input options\n";
print($OFILE join($OFS,("Condition","TsTv_A","TsTv_B","TsTv_filt","TsTv_filtN","TsTv_filtNAB")),"\n");

foreach my $cond (keys %fullconds)
{
    foreach my $fullcond (keys %{$fullconds{$cond}})
    {
        print($OFILE join($OFS,($fullcond,$Aconds{$cond},$Bconds{$cond},$filtconds{$cond},$filtNconds{$cond},${$fullconds{$cond}}{$fullcond})),"\n");
    }
}
close($OFILE);
exit;
