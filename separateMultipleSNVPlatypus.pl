use warnings;
use strict;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#Config
my $usage="Usage: $0 -i input_file -o output_file\nThis script reads a vcf file generated with platypus and splits multi-SNV lines in independent single-SNV lines\n";
my $description_VCF="##Multi-SNV variants separated in different entries using separateMultipleSNVPlatypus.pl";
my $info_header_line='##INFO=<ID=MultiSNV,Number=1,Type=String,Description="|-separated alternative alleles of the original multi-SNV VCF line">';
my $info_message="MultiSNV=";
my %forbidden_infos=map {$_=>1} qw(GL);

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
			
			##Detect fields of sample that should not be splitted based on the format (comma-separated but not related to multiSNVs)
			$ref_forbidden_formats=get_forbidden_format($columns[8],\%forbidden_infos);

			for ($i=0; $i<$nout; ++$i)
			{
				#Push constants and alt
				@out=();
				push(@out,@columns[0..3],$alts[$i],@columns[5..6]);
				push(@out,get_field_alt_addinfo(\@in_info,$i,";",\@alts));
				##Add format
				push(@out,$columns[8]);
				##Add samples
				for ($j=9; $j< scalar @columns; ++$j)
				{
					push(@out,get_field_alt_forbidden([split(":",$columns[$j],$j)],$i,":",$ref_forbidden_formats));
				}
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

sub get_field_alt_forbidden
{
	my ($ref_input,$i,$char,$ref_skip_i)=@_;
	my @columnout;
	my @altsthiscol;
	my $code;

	#print("DEBUG: @{$ref_input}\n");
	#print("DEBUG: @{$ref_skip_i}\n");
	for (my $j=0; $j<scalar @{$ref_input}; ++$j)
	{
	    if($ref_skip_i->[$j]==0 && $ref_input->[$j] =~m/,/)
	    {
	        @altsthiscol=split(",",$ref_input->[$j]);
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
	        push(@columnout,$ref_input->[$j]);
	    }
	}
	return join($char,@columnout);
}

sub get_forbidden_format
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
	my ($ref_input,$i,$char,$ref_alts)=@_;
	my $info_to_add=$info_message.join("|",@{$ref_alts});
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
