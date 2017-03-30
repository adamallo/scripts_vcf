use strict;
use warnings;

if (scalar @ARGV != 1 or $ARGV[0] eq "help" or ! -d $ARGV[0])
{
	print("Usage $0: directory\n");
	exit(1);
}

my @dirs = grep { -d } glob "$ARGV[0]/*";
my @samples=("A","B");
my @filters=("N","NAB");

for my $dir (@dirs)
{
	for my $sample (@samples)
	{
		for my $filter (@filters)
		{
			##Input. This is quite messy and hardcoded... Oh well, I am in a rush.
			my @filesvariants=glob("$dir/${sample}filt${filter}#*different.vcf.annotated.variant_function");
			scalar @filesvariants != 1 and die "Error filtering the files";
			my $filevariants=$filesvariants[0];
			my $filecounts="$dir/${sample}_${filter}_counts.csv";
			-f $filecounts or die "Error filtering the files. The file $filecounts does not exist";
			
			##Get annovar variants
			open(my $AVARIANTS,$filevariants);
			my @avariants=<$AVARIANTS>;
			close($AVARIANTS);
			
			##Dictionary of annotated variants	
			my %avardict=();
			open(my $OANNR,">$dir/${sample}_${filter}.INDELS");
			
			foreach my $var (@avariants)
			{
				next unless $var =~ /\S/;
				my @cols=split("\t",$var);
				my $id="$cols[2]_$cols[3]";
				if (length($cols[5]) != 1)
				{
					print($OANNR $var); ##We do not treat INDELS
				}
				else
				{
	#				print("DEBUG: $id\n");
					$avardict{$id}=$var;
				}
			}
	
			close($OANNR);
			
			open(my $COUNTS,$filecounts);
			my @counts=<$COUNTS>;
			close($COUNTS);
	
			##Dictionary of checked positions
			my %countsdict=();
			foreach my $var (@counts)
			{      
				next unless $var =~ /\S/;
                                next if $var =~/^#/;##Comment
				my @cols=split("\t",$var);
				print("WTF: $cols[4]\n");
	                        my $id="$cols[0]_$cols[1]";
				my $p="NA";
#				print("DEBUG: $id\n");
				if ($cols[4] != 0)
				{
					$p=$cols[5]/$cols[4];
				}
	                        $countsdict{$id}=[$cols[4],$cols[5],$p];
			}

			open(my $OANN,">$dir/${sample}_${filter}.SNVs");	
			foreach my $id (keys %avardict)
			{
				my $toadd=join("\t",$countsdict{$id});
				print($OANN "$avardict{$id}\t$toadd");
			}			
			close($OANN);
		}		
	}
}
