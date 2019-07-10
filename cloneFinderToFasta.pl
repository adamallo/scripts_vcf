use warnings;
use strict;

#Configuration variables
my $OFS="\t";
my $IFS="\t";
my $MFS=",";

#Flags
my $help;

my $usage="Usage: $0 inputfile.meg inputfile.tsv > outputfile\n The first input file corresponds to the *.meg file that clonefinder generates. The second input file is the original input file of clonefinder\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
if (scalar @ARGV < 2 || scalar @ARGV > 3 ||! -f $ARGV[0] || ! -f $ARGV[1])
{
	die "ERROR: problem parsing input files.\n$usage";
}

mopen(my $FH_INFILE,$ARGV[0]);
my @clonedata=<$FH_INFILE>;
close($FH_INFILE);

mopen($FH_INFILE,$ARGV[1]);
my @inputdata=<$FH_INFILE>;
close($FH_INFILE);

my $filterInvariants=0;
if (scalar @ARGV >2)
{
	$filterInvariants=$ARGV[2];
}

shift @inputdata;
my @vardata;
my @ref;
my @alt;
foreach my $line (@inputdata)
{
	chomp($line);
	@vardata=(split($IFS,$line))[1,2];
	push(@ref,$vardata[0]);
	push(@alt,$vardata[1]);
}
my $iclone=0;
my %outcontent;
my @clones;
foreach my $line (@clonedata)
{
	chomp($line);
	if ($line =~ /^[\#\! ]/)
	{
		next;
	}
	else
	{
		$outcontent{">Clone$iclone"}=translateVariants($line,\@ref,\@alt);
		push(@clones,">Clone$iclone");
		$iclone++;
	}
}

$outcontent{">Normal"}=[@ref];
push(@clones,">Normal");

my @mask;

if ($filterInvariants == 1)
{
	foreach my $clone (@clones)
	{
		for (my $pos=0; $pos< scalar @ref; ++$pos)
		{
			if( defined $mask[$pos] && $mask[$pos] ne $outcontent{$clone}->[$pos])
			{
				$mask[$pos]="0";
			}
			else
			{
				$mask[$pos]=$outcontent{$clone}->[$pos]
			}
		}
	}
}
else
{
	@mask=(0) x scalar @ref;
}
for (my $iclone=0; $iclone< scalar@clones; ++$iclone)
{
	print($clones[$iclone]."\n");
	for (my $pos=0; $pos< scalar @ref; ++$pos)
	{
		if($mask[$pos] eq "0")
		{
			print($outcontent{$clones[$iclone]}->[$pos]);
		}
	}
	print("\n");
}

exit;

sub translateVariants
{
	my @outline;
	my ($inline,$refs,$alts)=@_;
	my @snvs=split("",$inline);

	for (my $pos=0; $pos< scalar @$refs; ++$pos)
	{
		push(@outline,$snvs[$pos] eq "A"?$refs->[$pos]:$alts->[$pos]);
	}

	return \@outline;
}


#My subrutines
##############

#sub write_FASTA_file
#{
#    my ($filename,$refdata,$refvars)=@_; 
#    mopen(my $OUT_MSA, ">$filename");
#
#    my @int_v;
#    my @sortedvars=nat_i_sorter{@int_v=split($OFS,$_);$int_v[0],$int_v[1]} keys %$refvars;
#    my @refs;
#    my @alts;
#    my @finalvars;
#    my @thisdata;
#
#    #Making a list of references and alternates to be used to print later
#    #Also filtering out INDELS  
#    foreach my $var (@sortedvars)
#    {
#        @thisdata=split($OFS,$var);
#        if(length $thisdata[2] == length $thisdata[3]) #INDEL filter
#        {
#            push(@finalvars,$var);
#            push(@refs,$thisdata[2]);
#            push(@alts,$thisdata[3]);
#        }
#    }
#
#    #Fake normal
#    print($OUT_MSA ">Normal\n");
#    for (my $ivar=0; $ivar<scalar @finalvars; ++$ivar)
#    {
#        print($OUT_MSA $refs[$ivar]);
#    }
#    
#    for (my $isample=0; $isample<$nsamples; ++$isample)
#    {
#        print($OUT_MSA "\n>$samplename[$isample]\n");
#        for (my $ivar=0; $ivar<scalar @finalvars; ++$ivar)
#        {
#            if($refdata->[$isample]->{$finalvars[$ivar]}->[0]=~/1/)
#            {
#                print($OUT_MSA $alts[$ivar]);
#            }
#            elsif($refdata->[$isample]->{$finalvars[$ivar]}->[0]=~/\?/)
#            {
#                print($OUT_MSA "?");
#            }
#            else
#            {
#                print($OUT_MSA $refs[$ivar]);
#            }
#        }
#    }
#    
#    close($OUT_MSA);
#}

##THESE SUBRUTINES SHOULD PROBABLY BE PART OF MY OWN MODULE/LIBRARY

# Parse vcf returning a hash reference with keys in the form CHROM$OFSPOS$OFSREF$OFSALT and value = array_ref [genotype (coded as 0/0, 0/1, 1/0, or 1/1),$nref+$nalt, $nref, $nalt]
##########################################################################
sub parse_vcf
{
    my ($vcf1_file)=@_;
    open(my $VCF1,$vcf1_file) or die "The file $vcf1_file cannot be opened";
    my @vcf1=<$VCF1>;
    close($VCF1);
    my $flag=0;
    my $i;
    my %hash;
    my $key;
    my $value;

    for ($i=0;$i<scalar @vcf1;$i++)
    {
        unless($flag==0 and $vcf1[$i]=~/^#/)
        {
            if ($flag==0)
            {
                $flag=1;
            }
            chomp($vcf1[$i]);
            $key=$value=$vcf1[$i];
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            $value=[(split(":",(split($IFS,$value))[9]))[0,4,5]];
            $value->[3]=$value->[2];
            $value->[2]=$value->[1]-$value->[3];
            $hash{$key}=$value;
        }
    }
    return \%hash;
}

# Parses a TSV file (covN or covB), returning a reference to a hash with a key indicating each variant and the value an array reference to [number of total reads, number of reference reads, number of alternative reads]
# TODO: ATTENTION: We can have multiple-SNV here. I think this should not be like that. I should eliminate the multiple-SNVs in the VCF file. Right now this would change all the covB filter structure and therefore I am patching it up here. If I fix it in the other place, I should remove the separation of mutiple alternatives here
sub parse_tsv
{
    my ($file)=@_;
    mopen(my $FT, $file);
    my @content=<$FT>;
    close($FT);
    my %hash;
    my ($chr, $pos, $ref, $alt, $nref, $nalt);
    foreach my $line (@content)
    {
        chomp($line);
        my @biallelicSNVs=separateTriallelic($line);
        foreach my $snv (@biallelicSNVs)
        {
            ($chr, $pos, $ref, $alt, $nref, $nalt)=split($IFS,$snv);
            $hash{join($OFS,$chr,$pos,$ref,$alt)}=[$nref+$nalt, $nref, $nalt];
        }
    }
    return \%hash;
}

# TODO: ATTENTION: We can have multiple-SNV here. I think this should not be like that. I should eliminate the multiple-SNVs in the VCF file. Right now this would change all the covB filter structure and therefore I am patching it up here. If I fix it in the other place, I should remove the separation of mutiple alternatives here
sub separateTriallelic
{
    my ($line)=@_;
    my @columns=split($IFS,$line);
    my @alts=split($MFS,$columns[3]);
    my @nalts=split($MFS,$columns[5]);
    scalar @alts != scalar @nalts and die "ERROR: tsv input with different number of alternative variants and information on the number of alternative reads\n";
    my @outcontent;
    for (my $i=0; $i<scalar @alts; ++$i)
    {
        push(@outcontent,join($IFS,@columns[0..2],$alts[$i],$columns[4],$nalts[$i]));
    }
    
    return @outcontent;
}

# Writes a vcf with the variables contained in a hash selected from another VCF file
# ##################################################################################

sub write_variant_vcf
{
    my ($ref_hash,$filename,$vcf,$comment)=@_;
    open(my $OFILE, ">$filename") or die "ERROR: The file $filename could not be opened for writting";
    open(my $IFILE, "$vcf");
    my @icontent= <$IFILE>;
    close($IFILE);
    my $flag=0;
    my %hash=%{$ref_hash};
    my $key;

    #Copying the header and adding a new line with filtering info
    #Then adding the variants that are present in the hash.
    for(my $i=0;$i< scalar @icontent; ++$i)
    {
        if ($flag==0 and $icontent[$i]=~/^##/)
        {
            print($OFILE $icontent[$i]);
        }
        elsif($flag==0)
        {
            print($OFILE "$comment\n$icontent[$i]");
            $flag=1;
        }
        else
        {
            $key=$icontent[$i];
            chomp($key);
            $key=~s/^([^\t]+)\t([^\t]+)\t[^\t]\t([^\t]+)\t([^\t]+)\t.*/$1$OFS$2$OFS$3$OFS$4/;
            #print("DEBUG: Key $key\n");
            if(exists $hash{$key})
            {
                print($OFILE $icontent[$i]);
                delete($hash{$key});
            }
            if(scalar keys %hash == 0)
            {
                last;
            }

        }

    }
    close $OFILE;
}

# Parses the sample names from a VCF and returns an array with them

sub get_sample_names_vcf
{
    my ($vcfname)=@_;
    my @names;
    
    mopen(my $FT, $vcfname);
    my @content=<$FT>;
    close($FT);
    
    foreach my $line (@content)
    {
        if($line =~m/^#[^#]/)
        {
            chomp($line);
            @names=split($IFS,$line);
            splice(@names,0,9);
            last;
        }
    }
    return @names;
}

# Open with error message
###################################################################################
sub mopen
{
    open($_[0],$_[1]) or die "ERROR: impossible to open the file $_[1] in ".($_[1]=~m/^>/?"write":"read")."mode\n";
    return $_[0];
}

