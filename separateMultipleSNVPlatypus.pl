use warnings;
use strict;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#Config
my $usage="Usage: $0 -i input_file -o output_file\nThis script reads a vcf file generated with platypus and splits multi-SNV lines in independent single-SNV lines\nATTENTION: GL and PL are kept as in the original VCF file";
my $description_VCF="##Multi-SNV variants separated in different entries using separateMultipleSNVPlatypus.pl";
my $info_header_line='##INFO=<ID=MultiSNV,Number=2,Type=String,Description="|-separated alternative alleles of the original multi-SNV VCF line, |-separated original genotypes">';
my $info_message="MultiSNV=";
my %forbidden_infos=map {$_=>1} qw(GL PL);
my %withref_infos=map {$_=>1} qw(AD);
my %genotype_infos=map {$_=>1} qw(GT);

#GETOPT
my ($input_file,$output_file)=("","");
my $help;

(! GetOptions(
    'input_file|i=s' => \$input_file,
    'output_file|o=s' => \$output_file,
    'help|h' => \$help,
                )) or (($output_file eq "") || ($input_file eq "") || $help) and die $usage;

##Reading input file
open(my $INPUT_FH, $input_file) or die "ERROR opening the input file $input_file\n";
my @content=<$INPUT_FH>;
close($INPUT_FH);

open(my $OUTPUT, ">$output_file") or die "ERROR opening the output file $output_file\n";

print("Working on the file $input_file, to generate $output_file\n");

##Main loop
my $is_content=0;
my $isinfo=0;
my @columns;
my @alts;
my $nout;
my $i;
my $j;
my @out;
my @in_info;
my $code;
my $ref_forbidden_formats;
my $ref_withref_formats;
my $ref_genotype_formats;
my $ref_out_genotypes; 

foreach my $line (@content)
{
	if($is_content==1)
	{
		if($line =~ m/^[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+,[^\t]+\t/) ##Multi-SNV line
		{
			##ADDING THE MAGIC HERE
			chomp($line);
			@columns=split("\t",$line);
			@alts=split(",",$columns[4]);
			$nout=scalar @alts;
			@in_info=split(";",$columns[7]);
			
			##Make masks for fields that must be treated in different ways
			$ref_forbidden_formats=get_isformat($columns[8],\%forbidden_infos);
			$ref_withref_formats=get_isformat($columns[8],\%withref_infos);
            $ref_genotype_formats=get_isformat($columns[8],\%genotype_infos);

			for ($i=0; $i<$nout; ++$i)
			{
				#Push constants and alt
				@out=();
                $ref_out_genotypes=[];
				push(@out,@columns[0..3],$alts[$i],@columns[5..6]);
                ##Not adding INFO yet, I cannot use push anymore
				##Add format
				$out[8]=$columns[8];
				##Add samples
				for ($j=9; $j< scalar @columns; ++$j)
				{
					$out[$j]=get_processed_sample_fields([split(":",$columns[$j],$j)],$i,":",$ref_forbidden_formats,$ref_withref_formats,$ref_genotype_formats,$ref_out_genotypes);
				}
				##Add INFO, now with the original genotypes
				#WORKING HERE. I NEED TO GET THE GENOTYPES IN THE FOR LOOP (MODIFYING THE GET_PROCESSED SUBRUTINE AND THEN MODIFY GET_FIELD...ADDINFO TO ADD THE GENOTYPES. THE FINAL FORMAT SHOULD BE ALT1|ALT2|ALTN,GENOTYPESAMPLE1|GENOTYPESAMPLE2|GENOTYPESAMPLEN
                $out[7]=get_field_alt_addinfo(\@in_info,$i,";",\@alts,$ref_out_genotypes);
				print($OUTPUT join("\t",@out),"\n");
			}
		}
		else ##Normal line, just print
		{
			print($OUTPUT $line);
		}
	}
	elsif($is_content==0 && $line =~m/^#[^#]/) ##Last VCF header line, add a comment and print
	{
		print($OUTPUT join("\n",$description_VCF,$line));
		$is_content=1;
	}
	elsif($isinfo==1 && !($line=~m/^##INFO=/))
	{
		print($OUTPUT join("\n",$info_header_line,$line));
		$isinfo=0;
	}
	else ##VCF Header, just print it
	{
		if($line=~m/^##INFO=/ )	
		{
			$isinfo=1;
		}
		print($OUTPUT $line);
	}
}

close($OUTPUT);
exit;

sub get_field_alt
{
	my ($ref_input,$i,$char)=@_;
	my @columnout;
	my @altsthiscol;
	my $code;

	foreach my $helpi (@{$ref_input})
	{
	    if($helpi =~m/,/)
	    {
	        @altsthiscol=split(",",$helpi);
	        $code=$altsthiscol[0];
	        if($code=~s/^([^=]+=).*$/$1/ == 0)
	        {
	            $code="";
	        }
	        $altsthiscol[0]=~s/^[^=]+=//g;
	        push(@columnout,$code.$altsthiscol[$i]);
	    }
	    else
	    {
	        push(@columnout,$helpi);
	    }
	}
	return join($char,@columnout);
}

#Returns the sample fields for the $i genotype in this multi-SNV line
#It uses two references to mask arrays (binary arrays), $ref_skip_i, $ref_withref_i
#ref_skip_i[$i]==1 are not processed (i.e., are copied as they are for all $i)
#ref_withref_i[$i]==1 are processed taking into consideration that the first position
#refers to the REF allele, and the rest to ALT alleles.
#It also writes into an array ref the genotype of this i
sub get_processed_sample_fields
{
	my ($ref_input,$i,$char,$ref_skip_i,$ref_withref_i,$ref_genotype_i,$ref_out_genotypes)=@_;
	my @columnout;
	my @altsthiscol;
    my @genotype;
	my $code;

	#print("DEBUG: @{$ref_input}\n");
	#print("DEBUG: @{$ref_skip_i}\n");
	for (my $j=0; $j<scalar @{$ref_input}; ++$j)
	{
	    if($ref_skip_i->[$j]==0 && $ref_input->[$j] =~m/,/) ##If this should not be skipped and has an array of data
	    {
	        @altsthiscol=split(",",$ref_input->[$j]);
	        $code=$altsthiscol[0]; #Get the code (if present) and first element
            
            if($ref_withref_i->[$j]==1) ##The first element is the reference. We keep it and then add the needed alternative
            {
                shift(@altsthiscol);
                push(@columnout,join(",",$code,$altsthiscol[$i]));
            }
            else
            {
    	        if($code=~s/^([^=]+=).*$/$1/ == 0) ##We get the code (if present, otherwise clean it from previous iteractions)
    	        {
    	            $code="";
    	        }
    	        $altsthiscol[0]=~s/^[^=]+=//g; ##Eliminate the code
	            push(@columnout,$code.$altsthiscol[$i]);
            }
	    }
        elsif($ref_genotype_i->[$j] == 1)
        {
            push(@{$ref_out_genotypes},$ref_input->[$j]);
            @genotype=split("/",$ref_input->[$j]);
            if($genotype[0]==$i+1)
            {
                if($genotype[1]==$i+1)
                {
                    push(@columnout,"1/1");
                }
                else
                {
                    push(@columnout,"1/0");
                }
            }
            elsif($genotype[1]==$i+1)
            {
                push(@columnout,"0/1");
            }
            else
            {
                push(@columnout,"0/0");
            }
        }
	    else
	    {
	        push(@columnout,$ref_input->[$j]);
	    }
	}
	return join($char,@columnout);
}

sub get_isformat
{
	my @output;
	my ($data,$ref_forbidden_infos)=@_;
	my @infos=split(":",$data);

	for (my $i=0; $i< scalar @infos; ++$i)
	{
		$output[$i]=exists $ref_forbidden_infos->{$infos[$i]}?1:0;
	}
	return \@output;
}

sub get_field_alt_addinfo
{
	my ($ref_input,$i,$char,$ref_alts,$ref_out_genotypes)=@_;
	my $info_to_add=$info_message.join("|",@{$ref_alts}).",".join("|",@{$ref_out_genotypes});
	my @columnout;
	my @altsthiscol;
	my $code;
	my $info_added=0;

	foreach my $helpi (@{$ref_input})
	{
		if($info_added==0 && ($info_to_add lt $helpi))
		{
			push(@columnout,$info_to_add);
			$info_added=1;
		}
	    if($helpi =~m/,/)
	    {
	        @altsthiscol=split(",",$helpi);
	        $code=$altsthiscol[0];
	        if($code=~s/^([^=]+=).*$/$1/ == 0)
	        {
	            $code="";
	        }
	        $altsthiscol[0]=~s/^[^=]+=//g;
	        push(@columnout,$code.$altsthiscol[$i]);
	    }
	    else
	    {
	        push(@columnout,$helpi);
	    }
	}
	return join($char,@columnout);
}
