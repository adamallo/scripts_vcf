use strict;
use warnings;

##Config
my $FS=",";

my $inputfile;
my $outputfile;

open(my $INPUT,$inputfile);
open(my $OUTPUT,">$outputfile");
my @content=<$INPUT>;


for (my $i=1; $i<scalar(@content);++$i)
{
  my @cells=split($FS,$content[$i]);
}
