use strict;
use warnings;


##Config
our $FS=",";
our $OFS=",";
my $usage="$0 input_file output_file";

#Functions
sub array_to_string
{
  my $out="";
  foreach my $data (@_)
  {
    $out="$out$data$OFS";
  }
  chop($out);
  return $out;
}

#Main program

((scalar(@ARGV) == 2) && (-f $ARGV[0])) or die $usage;

my $inputfile=$ARGV[0];
my $outputfile=$ARGV[1];

#IO init
open(my $INPUT,$inputfile);
open(my $OUTPUT,">$outputfile");
my @content=<$INPUT>;
my @header=split($FS,$content[0]);
splice(@header,1,1);
print($OUTPUT array_to_string(@header));

#Main loop
for (my $i=2; $i<scalar(@content);$i+=2)
{
  my @cells1=split($FS,$content[$i-1]);
  my @cells2=split($FS,$content[$i]);
  my @outcells;
  scalar(@cells1) != scalar (@cells2) and die "Rows have different number of cells!";
  splice(@cells1,1,1);
  @outcells=splice(@cells1,0,2);
  splice(@cells2,0,3);
  for (my $j=0; $j<scalar(@cells1);++$j)
  {
    if ($cells1[$j] =~ /NA/ or $cells2[$j] =~ /NA/)
    {
      $outcells[$j+2]="NA";
    }
    else
    {
      print("DEBUG: $cells1[$j], $cells2[$j]\n");
      $outcells[$j+2]=abs($cells1[$j]-$cells2[$j]);
    }
  }
  print($OUTPUT array_to_string(@outcells),"\n");
}
exit;
