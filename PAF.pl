#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::HTS::Tabix;

my $OFS=",";
my $GNOMAD=$ENV{'GNOMAD'};
my $MAX_OPEN_IT=10;

my $usage="\n$0 variants.vcf PAF threshold.\n-------------------------------\nThis program calculates the proportion of variants with a population allele frequency less or equal than a user-specified threshold. It requires the GNOMAD summary file to be indicated by the environment variable \$GNOMAD\n\n";

#ARGV
if(scalar @ARGV != 2 || ! -f $ARGV[0] || !($ARGV[1]>=0 && $ARGV[1]<=1))
{
    print $usage;
    exit 1;
}

my $ifile=$ARGV[0];
my $paf=$ARGV[1];

print(STDERR "Parsing $ifile...");

my $ref_variants=parse_vcf($ifile);

print(STDERR "\tDone\n");

print(STDERR "Calculating the percentage of variants in $ifile that have a PAF <= $paf...");

my $ref_pafdata=getPAFdata($ref_variants);

my $tpaf=0;
my $total=scalar keys %$ref_pafdata;
my $pass=0;
foreach my $var (keys %$ref_pafdata)
{
    $tpaf=$ref_pafdata->{$var}[1];
    if ($tpaf eq "NA" || $tpaf<= $paf)
    {
        $pass+=1;
    }
}
print(STDERR "\tDone\n");

if ($total>0)
{
    print(($pass*1.0)/$total);
}
else
{
    print("NaN");
}

print("\n");
exit 0;

sub getPAFdata
{
    my $tabix = Bio::DB::HTS::Tabix->new( filename =>$GNOMAD );
    my %outdata; #key: CHROM${OFS}POS${OFS}REF${OFS}ALT value: [$tfilt,$taf];
    my $parsedPAFS;
    my %obtainedPAFs; #key = CHROM${OFS}POS${OFS}REF${OFS}ALT #value = [$tfilt, $taf]
    my $tabix_iter;
    my ($chr, $nstart, $ref, $alt);
    my $line;
    my $i;
    my ($tstart,$talt,$tref,$tfilt,$taf);

    foreach my $ref_hash (@_)
    {
        foreach my $key (keys %{$ref_hash}) ##Foreach key
        {
            #key = CHROM${OFS}POS${OFS}REF${OFS}ALT
            unless (exists $obtainedPAFs{$key})
            {
                ##Obtain gnomAD data for this genomic position if it has not been obtained before. MULI-SNVs are represented in different lines
                ($chr,$nstart,$ref,$alt)=split($OFS,$key);
                $tabix_iter=$tabix->query("$chr:$nstart-".($nstart+1));
                if(defined $tabix_iter)
                {
                    while($line=$tabix_iter->next)
                    {
                        #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  S1 ... SN
                        $line =~ s/^[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t[^\t]*AF=([^;]+).*$/$1\t$3\t$2\t$4\t$5/;
                        ($tstart,$talt,$tref,$tfilt,$taf)=split("\t",$line);
                        $obtainedPAFs{"$chr$OFS$tstart$OFS$tref$OFS$talt"}=[$tfilt,$taf];
                    }

                }
#DEBUG
#                else
#                {
#                    warn "There is no data for $chr:$nstart-".($nstart+1);
#                }
            }

            unless (exists $outdata{$key}) ##If it exists we don't need to do anything, otherwise, we copy the data if we have it, or add NA data since we know we have tried to obtain it and it is not available
            {
                if(exists $obtainedPAFs{$key})
                {
                    $outdata{$key}=$obtainedPAFs{$key};
                }
                else
                {
                    $outdata{$key}=["NA","NA"];
                }
            }
        }##foreach key
    }##foreach hash
    return \%outdata;
}

# Parse vcf, splitting multi-snv lines in different lines
#########################################################
sub parse_vcf
{
    my ($vcf1_file)=@_;
    my $VCF1;
    myopen(\$VCF1,$vcf1_file) or die "The file $vcf1_file is not located in the output folder. Please, check if the variant caller has successfully finished";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            ##KEY will be CHROM$OFSPOS$OFSREF$OFSALT
            $key=$vcf1[$i];
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            chomp($key);
            $hash{$key}=1;
#print("DEBUG: New variant being hashed $key\n");
        }
    }
    return \%hash;
}

#Retries the opening of a file a number of times, instead of failing at the first attempt
#This seems to be needed when running on a network attached system in parallel with many processes working on it
sub myopen
{
    my ($refFH,$name)=@_;
    my $i=0;
    while (!open($$refFH,$name))
    {
        warn "Retrying to open the file $name, attempt $i\n";
        sleep 10;
        $i+=1;
        if ($i > $MAX_OPEN_IT)
        {
            return 0;
        }
    }
    return 1;
}
