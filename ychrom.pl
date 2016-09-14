use strict;
use warnings;

use Bio::DB::Sam;


##IO
my $usage="$0 bamfile FASTA(refgenome)\n";
my $bamfile="";
my $fasta="";

if (scalar @ARGV != 2 || ! -f $ARGV[0] || ! -f $ARGV[1])
{
  die $usage;
}

##Conf
my $mapq=20;
my $phred=20;


my $sam = Bio::DB::Sam->new(-bam  => $bamfile, -fasta=> $fasta);
my @alignments = $sam->get_features_by_location(-type=>'read_pair', -seq_id => 'chrY');
my @out_alignments;
foreach my $alignment (@alignments)
{
  my $filter=0;
  if($alignment->qual >= 5) #MAPQ
  {
    #Phred scores, commented out for now
    # my @baseqs=$alignment->qscore
    # foreach my $baseq (@baseqs)
    # {
    #   if ($baseq<$phred)
    #   {
    #     $filter=1;
    #     last;
    #   }
    # }
    unless ($filter)
    {
      push(@out_alignments,$alignment);
    }
  }

}

print(scalar @out_alignments," pared reads aligned to the Y chromosome with MAPQ >= $mapq\n");
