#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

my $outfile="";
our $OFS="\t";
my $usage="$0 main_directory output_file\n";

if (scalar @ARGV != 2 || !-d $ARGV[0])
{
    die($usage);
}

open(my $OUTPUT,">$ARGV[1]") or die "Impossible to open the output file $ARGV[1]\n";

print($OUTPUT "Patient${OFS}Common${OFS}Common_pos${OFS}Common_type${OFS}A${OFS}A_pos${OFS}A_type${OFS}B${OFS}B_pos${OFS}B_type\n");

chdir($ARGV[0]);

opendir(my $DH,".");

my @dirs = grep {-d "\./$_" && ! /^\.{1,2}$/} readdir($DH);

closedir($DH);

foreach my $dir (@dirs)
{
    my $name=basename($dir);
    my @common_files=glob("$dir/filtNABU#*common.vcf.annotated.variant_function");
    my $common_file=$common_files[0];
    my $common_file_aux=$common_file;
    $common_file_aux=~s/variant_function/exonic_variant_function/;

    if (scalar @common_files != 1)
    {
#        print("DEBUG: @common_files\n");
        die "More than one files detected with the format filtNABU#*common.vcf.annotated.variant_function. This script is intended to analyse the output of only one filter\n";
    }

    my @da_files=glob("$dir/AfiltNAB#*different.vcf.annotated.variant_function");
    my $da_file=$da_files[0];
    my $da_file_aux=$da_file;
    $da_file_aux=~s/variant_function/exonic_variant_function/;

    if (scalar @da_files != 1)
    {
        die "More than one files detected with the format AfiltNAB#*different.vcf.annotated.variant_function. This script is intended to analyse the output of only one filter\n";
    }

    my @db_files=glob("$dir/BfiltNAB#*different.vcf.annotated.variant_function");
    my $db_file=$db_files[0];
    my $db_file_aux=$db_file;
    $db_file_aux=~s/variant_function/exonic_variant_function/;

    if (scalar @db_files != 1)
    {
        die "More than one files detected with the format BfiltNAB#*different.vcf.annotated.variant_function. This script is intended to analyse the output of only one filter\n";
    }

    open(my $COMMON,$common_file) or die "Impossible to open the input file $common_file\n";
    open(my $COMMON_AUX,$common_file_aux) or die "Impossible to open the input file $common_file_aux\n";
    open(my $DA, $da_file) or die "Impossible to open the input file $da_file\n";
    open(my $DA_AUX, $da_file_aux) or die "Impossible to open the input file $da_file_aux\n";
    open(my $DB, $db_file) or die "Impossible to open the input file $db_file\n";
    open(my $DB_AUX, $db_file_aux) or die "Impossible to open the input file $db_file_aux\n";

    my @c_content=<$COMMON>;
    my @c_content_aux=<$COMMON_AUX>;
    my @da_content=<$DA>;
    my @da_content_aux=<$DA_AUX>;
    my @db_content=<$DB>;
    my @db_content_aux=<$DB_AUX>;

    my @lengths=(scalar @c_content,scalar @da_content,scalar @db_content);
    @lengths=sort(@lengths);
   
    close($COMMON);
    close($DA);
    close($DB);
    close($COMMON_AUX);
    close($DA_AUX);
    close($DB_AUX);
    my @temp_data;
    my $type="";
    my @temp2=();

    for (my $i=0; $i<$lengths[2];++$i)
    {
        my @out_line=("$name");
        if(scalar(@c_content) <= $i)
        {
            @temp_data=("","","");
        }
        else
        {
            @temp_data=split("\t",$c_content[$i]);
            if($temp_data[0]=~m/^exonic$/)
            {
                foreach my $line (@c_content_aux)
                {
                    my $n=$i+1;
                    if($line=~/^line$n/)
                    {
                        @temp2=split("\t",$line);
                        $type=$temp2[1];
                        last;
                    }
                }
            }
            else
            {
                $type="";
            }
            @temp_data=($temp_data[1],$temp_data[0],$type);
        }
        #print("DEBUG: @temp_data\n");
        push(@out_line,@temp_data);
        #print("DEBUG: @out_line\n");
        if(scalar(@da_content) <= $i)
        {
            @temp_data=("","","");
        }
        else
        {
            @temp_data=split("\t",$da_content[$i]);
           if($temp_data[0]=~m/^exonic$/)
            {
                foreach my $line (@da_content_aux)
                {
                    my $n=$i+1;
                    if($line=~/^line$n/)
                    {
                        @temp2=split("\t",$line);
                        $type=$temp2[1];
                        last;
                    }
                }
            }
            else
            {
                $type="";
            }

             @temp_data=($temp_data[1],$temp_data[0],$type);
        }
        push(@out_line,@temp_data);
        
        if(scalar(@db_content) <= $i)
        {
            @temp_data=("","","");
        }
        else
        {
            @temp_data=split("\t",$db_content[$i]);
            if($temp_data[0]=~m/^exonic$/)
            {
                foreach my $line (@db_content_aux)
                {
                    my $n=$i+1;
                    if($line=~/^line$n/)
                    {
                        @temp2=split("\t",$line);
                        $type=$temp2[1];
                        last;
                    }
                }
            }
            else
            {
                $type="";
            }

            @temp_data=($temp_data[1],$temp_data[0],$type);
        }
        #print("DEBUG: @temp_data\n");
        push(@out_line,@temp_data);
#        print("DEBUG: @out_line\n");
        print($OUTPUT print_array(@out_line),"\n");
    }
}

close($OUTPUT);
exit;

sub print_array
{
    my @array=@_;
    my $output="";
    foreach my $element (@array)
    {
        $output="$output$element$OFS";
    }
    chop($output);
    return($output);
}
