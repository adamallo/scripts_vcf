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

#Hash of hashes keys: common filtering conditions, value hash
#   secondary hashes: key=full filtering condition, value = tstv
my %fullUconds;
my %fullIconds;

#Hashes, keys=filtering_condition, value=tstv
my %Aconds;
my %Bconds;
my %filtUconds;
my %filtIconds;
my %filtNUconds;
my %filtNIconds;

my @temp;
my $cond;

foreach my $line (@content)
{
    chomp($line);
    @temp=split($FS,$line);
    if($temp[0]=~m/different/)
    {
        next;
    }
    if($temp[0]=~s/^filtNABU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $cond=(split("${sep_param}NAB${sep_param}",$temp[0]))[0];
        if(defined $fullUconds{$cond})
        {
            #print("DEBUG $cond, $fullconds{$cond}\n");
            ${$fullUconds{$cond}}{$temp[0]}=$temp[1];
        }
        else
        {
            my %hash;
            $hash{$temp[0]}=$temp[1];
            $fullUconds{$cond}=\%hash;
        }
    }
    elsif($temp[0]=~s/^filtNABI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $cond=(split("${sep_param}NAB${sep_param}",$temp[0]))[0];
        if(defined $fullIconds{$cond})
        {
            #print("DEBUG $cond, $fullconds{$cond}\n");
            ${$fullIconds{$cond}}{$temp[0]}=$temp[1];
        }
        else
        {
            my %hash;
            $hash{$temp[0]}=$temp[1];
            $fullIconds{$cond}=\%hash;
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
    elsif($temp[0]=~s/^filtNU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtNUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtNI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtNIconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtIconds{$temp[0]}=$temp[1];
    }

}

open(my $OFILE, ">$outputfile") or die "The file $outputfile can't be opened, check your input options\n";
print($OFILE join($OFS,("Condition","TsTv_A","TsTv_B","TsTv_filtU","TsTv_filtI","TsTv_filtNU","TsTv_filtNI","TsTv_filtNABU","TsTv_filtNABI")),"\n");

###Pending: I should make sure that I and U have the same keys for every filtering level

foreach my $cond (keys %fullUconds)
{
    foreach my $fullcond (keys %{$fullUconds{$cond}})
    {
        print($OFILE join($OFS,($fullcond,$Aconds{$cond},$Bconds{$cond},$filtUconds{$cond},$filtIconds{$cond},$filtNUconds{$cond},$filtNIconds{$cond},${$fullUconds{$cond}}{$fullcond},${$fullIconds{$cond}}{$fullcond})),"\n");
    }
}
close($OFILE);
exit;
