use strict;
use warnings;

use Bio::DB::Sam;
use Bio::DB::Sam::Constants;
use Bio::FeatureIO;
use Data::Dumper;

##Conf
our $mapq=20;
our $phred=20;

#Data
my $par1_start=10001;
my $par1_end=2649520;
my $par3_start=2926383;
my $par3_end=6543100;
my $par2_start=59034050;
my $par2_end=59363566;

#Subrutines
###########
sub filter
{
  my $read=shift;
  my $filter=0;
  if($read->get_tag_values(FLAGS->{0x002}) == 1)
  {
    if($read->qual >= $mapq) #MAPQ
    {
      #Phred scores, commented out for now
#       my @baseqs=$read->qscore;
#       foreach my $baseq (@baseqs)
#       {
#         if ($baseq<$phred)
#         {
#           $filter=0;
#           last;
#         }
#       }
      $filter=1;
    }
  }
  return $filter;

}

sub substract_regions
{   
    my ($ref_final_regions,$ref_excluded_regions)=@_;
    my @excluded_regions=@{$ref_excluded_regions};
    my @final_regions=@{$ref_final_regions};
    my $start_e=0;
    my $end_e=0;
    my $start_o=0;
    my $end_o=0;
    
    foreach my $reg_to_exclude (@excluded_regions)
    {
        $start_e=${$reg_to_exclude}[0];
        $end_e=${$reg_to_exclude}[1];
#        print("DEBUG: reg to exclude $start_e, $end_e\n");
        my @original_final_regions=@final_regions;
        @final_regions=();
        for (my $i=0;$i<scalar(@original_final_regions);++$i)
        {
#            print("\tExternal loop, original regions=",scalar @original_final_regions,"\n");
            #my $reg=$final_regions[$i];
            my $reg=$original_final_regions[$i];
            $start_o=${$reg}[0];
            $end_o=${$reg}[1];
#            print("\tDEBUG: modifying region $start_o, $end_o\n");
            if($start_e<=$start_o && $end_e>=$end_o) ###eliminate this element of the array
            {
#                print("\t\tDEBUG: removing this original region (not copying it)");
                #splice(@final_regions,$i,1);
            }
            elsif ($start_e>$start_o && $end_e<$end_o) ###Split segment in two
            {
#                print("\t\tDEBUG: Split old: ",${$reg}[0],",",${$reg}[1],"; new: ",${$reg}[0],",",$start_e-1,";","$end_e+1, $end_o\n");
                my @new=($end_e+1,$end_o);
                ${$reg}[1]=$start_e-1;
                push(@final_regions,$reg);
                push(@final_regions,\@new);
                #splice(@final_regions,$i+1,0,\@new);
            }
            elsif ($start_e<=$start_o && $end_e>$start_o)##Prune from the left
            {
#                print("\t\tDEBUG: changing the left side of the fragment from ",${$reg}[0],"to $end_e+1\n");
                ${$reg}[0]=$end_e+1;
                push(@final_regions,$reg);
            }
            elsif($end_e>=$end_o && $start_e<$end_o)##Prune from the right
            {
#                print("\t\tDEBUG: changing the right side of the fragment from",${$reg}[1],"to $start_e-1\n");
                ${$reg}[1]=$start_e-1;
                push(@final_regions,$reg);
            }
            else
            {
                push(@final_regions,$reg);
            }
#            print("\tDEBUG: final regions= ",scalar @final_regions,"\n");
        }
    }
    
    #DEBUG
#    foreach my $reg (@final_regions)
#    {
#        print("start= ",${$reg}[0],"end= ",${$reg}[1],";");
#    }
#    print("\n");
    return(@final_regions);
}

##IO
my $usage="$0 bamfile FASTA(refgenome) bed_file min_mapq\n\n
This script calculates the percentage of the total well aligned read pairs (0x008 SAM tag) over certain mapping quality threshold that map to the Y chromosome\n";
my $bamfile="";
my $fasta="";
my $bedfile="";

if (scalar @ARGV != 4 || ! -f $ARGV[0] || ! -f $ARGV[1] || ! -f $ARGV[2])
{
  die $usage;
}
else
{
    $bamfile=$ARGV[0];
    $fasta=$ARGV[1];
    $bedfile=$ARGV[2];
    $mapq=$ARGV[3];
}

my $sam = Bio::DB::Sam->new(-bam  => $bamfile, -fasta=> $fasta, -expand_flags=>1);
my @ids = $sam->seq_ids();

#my @total = $sam->get_features_by_location(-type=> 'read_pair');
#my $unfiltered_total_pairs=scalar(@total);
#@total=();
#print("Total, memory usage should drop now\n");


##Total filtered reads
######################
my $filtered_total_pairs=0;

foreach my $region (@ids)
{
    my @filtered_total = $sam->get_features_by_location(-type=>'read_pair',-seq_id=>$region,-filter => \&filter);
    $filtered_total_pairs+=scalar(@filtered_total);
#    print("DEBUG: adding ",scalar(@filtered_total),"alignments to the total\n");
    @filtered_total=();
}


#my @pairs = $sam->get_features_by_location(-type=>'read_pair',-seq_id => 'Y');
#my $unfiltered_y_pairs=scalar(@pairs);
#@pairs=();

##Filtered reads chromosome Y
#############################
my @filtered_y = $sam->get_features_by_location(-type=>'read_pair',-seq_id => 'Y',-filter => \&filter);
my $filtered_y_pairs=scalar(@filtered_y);
#@filtered_y=();

##Filtered reads non-PAR chromosome Y
#####################################

#Defining the regions to look for variants
my @excluded_regions=([$par1_start,$par1_end],[$par3_start,$par3_end],[$par2_start,$par2_end]); #ordered array of anonymous arrays (array refs generated by [])
my @final_regions=([1,$sam->length("Y")]);
@final_regions=substract_regions(\@final_regions,\@excluded_regions);

my @filtered_nonpar;
foreach my $region (@final_regions)
{
    my $start=${$region}[0];
    my $end=${$region}[1];
    push(@filtered_nonpar,$sam->get_features_by_location(-type=>'read_pair',-seq_id => 'Y',-start=>$start,-end=>$end,-filter => \&filter));
    #print("DEBUG: adding another region of filtered nonpar, current number of alignments ",scalar(@filtered_nonpar),"\n");
}
my $filtered_nonpar_pairs=scalar @filtered_nonpar;
#print("$unfiltered_y_pairs read pairs aligned to the Y chromosome\n");


##Filtered reads non-PAR BED chromosome Y
#########################################

my $bed = Bio::FeatureIO->new(-file => $bedfile, -format => 'bed');
my @bed_regions;

#print("Debug:entering the loop\n");

while (my $entry = $bed->next_feature())
{
    if ($entry->seq_id =~ m/Y/)
    {
        push(@bed_regions,[$entry->start,$entry->end]);
    }
}

#print("DEBUG: BED regions: ");
#foreach my $ref (@bed_regions)
#{
#    print(${$ref}[0],", ",${$ref}[1],"; ");
#}
#print("\n");

@bed_regions=substract_regions(\@bed_regions,\@excluded_regions);

my @filtered_bed_nonpar;
foreach my $region (@bed_regions)
{
    my $start=${$region}[0];
    my $end=${$region}[1];
    push(@filtered_bed_nonpar,$sam->get_features_by_location(-type=>'read_pair',-seq_id => 'Y',-start=>$start,-end=>$end,-filter => \&filter));
    #print("DEBUG: adding another region of filtered nonpar, current number of alignments ",scalar(@filtered_nonpar),"\n");
}
my $filtered_bed_nonpar_pairs=scalar @filtered_bed_nonpar;

#print("Total,CHRY,CHRY-NONPAR,CHRY-NONPAR-BED\n");
print("$filtered_total_pairs,$filtered_y_pairs,$filtered_nonpar_pairs,$filtered_bed_nonpar_pairs\n");
#print("$filtered_y_pairs paired reads aligned to the Y chromosome with MAPQ >= $mapq\n");
#print($filtered_y_pairs/($filtered_total_pairs+0.0)," percentage of total reads with MAPQ >= $mapq in Y \n");
