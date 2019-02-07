#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;
use Env;
use Bio::DB::HTS::Tabix;

my $usage="Usage: $0 [options] -i input_file -o output_file -f main_folder [-s suffix]\nThe variants file is the 1col result of tabulate_annovar_angelo.pl, filtered for the desired variants\nThis script adds the total number of reads, the number of alternative reads and the population allele frequency for each variant\nThe vcf parsed are in the form of A|B\${suffix}.*[different|common]. It requires the environment variable GNOMAD pointing to the GNOMAD database\n";
######################################################

######################################################
##MAIN
######################################################
##Configuration variables
my $IFS="\t";
my $OFS="\t";
my $FFS=","; #Separator of two cells of the same content

##Input variables
my $inputDir="";
my $output_file="";
my $input_file="";
my $suffix="filtcovBNABPAF";
my $help;

##Getopt
######################################################
(! GetOptions(
    'input_file|i=s' => \$input_file,
    'output_file|o=s' => \$output_file,
    'main_folder|f=s' => \$inputDir,
    'suffix|s=s' =>\$suffix ,
    'help|h' => \$help,
                )) or (($output_file eq "") || ($inputDir eq "")  || (! -s $input_file) || $help) and die $usage;


#Parse input_file
open(my $IFILE, $input_file);
my @variantsContent=<$IFILE>;
close($IFILE);

#Get environment variables
my $GNOMAD=$ENV{'GNOMAD'};

my @row;
my $case;
my $shared;
my $chr;
my $start;
my $end;
my $ref;
my $alt;
my $varid;
my $gene;
my $kind;
my $type;
my %variantdata;
my $header=$variantsContent[0];

print("Parsing input tabulated data...");
for (my $i=1; $i<scalar @variantsContent; ++$i){ ##Skips the header
    @row=split("$IFS",$variantsContent[$i]);
    scalar @row != 10 and die "The line $i of the variant file, $input_file does not conform with the expected format since it has ".scalar @row." columns instead of 10\n";
    ($case,$shared,$chr,$start,$end,$ref,$alt,$gene,$kind,$type)=@row;
    chomp($row[9]);
    $varid="${chr}$OFS${start}$OFS${ref}$OFS${alt}";

    if (exists $variantdata{$case}){
        if(exists $variantdata{$case}{$shared}){
            if(exists $variantdata{$case}{$shared}{$varid}){
                die "Error, several entries for the same variant";
            }else{
                $variantdata{$case}{$shared}{$varid}=join($OFS,@row);
            }
        }
        else{
            $variantdata{$case}{$shared}={$varid => join($OFS,@row)};
        }
    }else{
        $variantdata{$case}={$shared => {$varid => join($OFS,@row)}};
    }

}
print(" Done\n");

print("Obtaining population allele frequency data for parsed variants...");
my $ref_PAF_data=getPAFdata(\%variantdata);
print(" Done\n");

##Main loop, directories
opendir(my $IDIR, $inputDir) or die "can't opendir $inputDir: $!";
my @dirs = grep { /^[^.]/ && -d "$inputDir/$_" } readdir($IDIR);
closedir $IDIR;

print("Writting the results...");
#my @shared_states=("Common","A","B");
##Opening the output file and prining the header
open(my $OUTPUT, ">$output_file") or die "Impossible to write the output file\n";
chomp($header);
$header=~s/$IFS/$OFS/g;
print($OUTPUT join("$OFS",$header,"VariantReads","TotalReads", "PAF"),"\n");

foreach my $dir (@dirs){

    $case=$dir;
    #$case=~s/^.*?_//g;
    #$case=~s///g;
    my @files=get_files($dir);   

    my $ref_a_private=parseVCF($files[0]);
    my $ref_b_private=parseVCF($files[1]);
    my $ref_a_common=parseVCF($files[2]);
    my $ref_b_common=parseVCF($files[3]);
    
    #print("DEBUG: ","case:$case",keys%{$variantdata{$case}->{"A"}});
    foreach my $varid (keys %{$variantdata{$case}->{"A"}})
    {
        exists $ref_a_private->{$varid} or die "Private A variant $varid not found in the private VCF file for case $case\n";
		print($OUTPUT join($OFS,$variantdata{$case}->{"A"}->{$varid},@{$ref_a_private->{$varid}},$ref_PAF_data->{$varid}[1]),"\n");
    }
    foreach my $varid (keys %{$variantdata{$case}->{"B"}})
    {
		exists $ref_b_private->{$varid} or die "Private B variant $varid not found in the private VCF file for case $case\n";
		print($OUTPUT join($OFS,$variantdata{$case}->{"B"}->{$varid},@{$ref_b_private->{$varid}},$ref_PAF_data->{$varid}[1]),"\n");
    }
    foreach my $varid (keys %{$variantdata{$case}->{"Common"}})
    {
		(exists $ref_a_common->{$varid} and exists $ref_b_common->{$varid}) or die "Common variant $varid not found in one or both common variant VCF files for case $case\n";
		print($OUTPUT join($OFS,$variantdata{$case}->{"Common"}->{$varid},join($FFS,$ref_a_common->{$varid}->[0],$ref_b_common->{$varid}->[0]),join($FFS,$ref_a_common->{$varid}->[1],$ref_b_common->{$varid}->[1]),$ref_PAF_data->{$varid}[1]),"\n");
    }
}

close($OUTPUT);
print("Done\n");
exit;

sub get_files {
    my $directory=$_[0];
    open(my $FILE, "$inputDir/$directory/vcfdict.csv");
    my @vcfdict_cont=<$FILE>;
    close($FILE);
    my @outfiles;
    my $detected=0;
    #print("DEBUG: ", join(" ",@vcfdict_cont));
    foreach my $file(@vcfdict_cont)
    {
        chomp($file);
        if($file =~ m/A$suffix.*different\.vcf/)
        {
            $outfiles[0]=$inputDir."/".$directory."/".(split(",",$file))[1];
            $detected+=1;
        }
        elsif($file =~ m/B$suffix.*different\.vcf/)
        {
            $outfiles[1]=$inputDir."/".$directory."/".(split(",",$file))[1];
            $detected+=1;
        }
    }

    ##The common files in the other output may not have the variant (one or the other). To get to the info we need to get to the ultimate source, the original vcf file. This is like this due to the fact that we are using the union of the variants.
    opendir(my $DIR, $inputDir."/".$directory) or die "can't opendir $directory: $!";
    my @initial_variant_files = grep { m/^[A|B](#--[^#]+)+.vcf$/ } readdir($DIR);
    closedir($DIR);
    if($initial_variant_files[0]=~m/^A/)
    {
        $outfiles[2]=$inputDir."/".$directory."/".$initial_variant_files[0];
        $outfiles[3]=$inputDir."/".$directory."/".$initial_variant_files[1];
        $detected+=2;
    }
    else
    {
        $outfiles[2]=$inputDir."/".$directory."/".$initial_variant_files[1];
        $outfiles[3]=$inputDir."/".$directory."/".$initial_variant_files[0];
        $detected+=2;
    }

    $detected != 4 and die "Problem parsing the VCF files. Detected $detected files\n";
    #print("DEBUG: @outfiles\n");
    return @outfiles;

}

sub parseVCF
{
    my $file=$_[0];
    open (my $FILE, $file) or die "Error opening the file $file\n";
    #print("DEBUG: file $file\n");
    my @rawdata=<$FILE>;
    close($file);
    my %data;
    my @temp;
    my @format;
    my @sample;
    my ($chr, $pos, $ref, $alt);
    my ($newalt,$newref,$newpos); #IMPORTANT NOTE: Annovar uses a different ref/alt format for INDELS. They never contain repeated information. For example, if ref is A and alt AC, in annovar this will be noted as - C.     ###MultiSNVS sometimes have extra common characters. For example, a reference AC, with alternatives CC and AT. When the multisnv is split in two variants, AC AT has an extra As, and annovar will consider it as C T. These two things generate downstream problems. I am adding a second entry with this format to solve it.
    my (@aref, @aalt);
    my $i;
    my $nshared;
    my $minlength;
    my ($nreads,$nvarreads);

    foreach my $line (@rawdata)
    {
        chomp($line);
        if (!($line =~ m/^#/))
        {
            @temp=split("\t", $line);
            ($chr, $pos, $ref, $alt)=@temp[0,1,3,4];
            @format=split(":", $temp[8]);
            @sample=split(":", $temp[9]);
            ($format[(scalar @format) -1] ne "NV") or ($format[(scalar @format) -2] ne "NR") and die "Format problem\n";
            ($nvarreads,$nreads)=@sample[scalar @sample -1, scalar @sample -2];
            $data{"$chr$OFS$pos$OFS$ref$OFS$alt"}=[$nvarreads,$nreads];

            @aref=split("",$ref);
            @aalt=split("",$alt);
            $nshared=0;
            $minlength=scalar @aref<scalar @aalt ? scalar @aref : scalar @aalt;
            
            for ($i=0; $i< $minlength; ++$i)
            {
                if($aref[$i] ne $aalt[$i])
                {
                    last;
                }
                else
                {
                    $nshared+=1;
                }
            }    
            
            if($nshared>0)
            {
                $newref=join("",splice(@aref,$nshared));
                $newalt=join("",splice(@aalt,$nshared));
                if($newref eq "")
                {
                    $newref="-";
                    $nshared-=1;
                }
                if($newalt eq "")
                {
                    $newalt="-";
                }
                
                $newpos=$pos+$nshared;
    
                $data{"$chr$OFS$newpos$OFS$newref$OFS$newalt"}=[$nvarreads,$nreads];
                #print("DEBUG: $chr$OFS$pos$OFS$ref$OFS$alt converted into $chr$OFS$newpos$OFS$newref$OFS$newalt\n");
            }
            

#            if($alt=~m/^$ref/) ##Adding a second entry, as explained in the "important note" right above
#            {
#                ($newref,$newalt)=($ref,$alt);
#                $newref="-";
#                $newalt=~s/^$ref//;
#                $data{"$chr$OFS$pos$OFS$newref$OFS$newalt"}=[$nvarreads,$nreads];
#            }
#            if($ref=~m/^$alt/) ##Adding a second entry, as explained in the "important note" right above
#            {
#                ($newref,$newalt)=($ref,$alt);
#                $newalt="-";
#                $newref=~s/^$alt//;
#                $newpos=$pos+length($alt);
#                #print("DEBUG: before $ref, $alt, $pos. after: $newref, $newalt, $newpos\n");
#                $data{"$chr$OFS$newpos$OFS$newref$OFS$newalt"}=[$nvarreads,$nreads];
#            }
            #print("DEBUG parseVCF: $chr$OFS$pos$OFS$ref$OFS$alt: $nvarreads,$nreads\n");
        }
    }
    #print("DEBUG:",join(",",keys %data),"\n");   
    return \%data;
}
 
##Returns a hash ref to the population allele frequency information for all variants present in the input reference hashes 
#The tabix query gets the data we are looking for and also some neighbors. To keep the output data clean from those, we handle two different hashes.
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
    my ($ref_hash)=@_;
    my ($tstart,$talt,$tref,$tfilt,$taf);

    foreach my $case (keys %{$ref_hash})
    {
        foreach my $shared (keys %{$ref_hash->{$case}})
        {
            foreach my $key (keys %{$ref_hash->{$case}->{$shared}})
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

            }#foreach varid
        }#foreach shared
    }#foreach case

    return \%outdata;
}

sub simplify_varid
{
    my ($id)=@_;
    my ($chr,$pos,$ref,$alt)=split($OFS,$varid);
    my @aref=split("",$ref);
    my @aalt=split("",$alt);
    my $nshared=0;
    
    for (my $i=0; $i< scalar @aref; ++$i)
    {
        if($aref[$i] ne $aalt[$i])
        {
            last;
        }
        else
        {
            $nshared+=1;
        }
    }    
    
    my $newref;
    my $newalt;
    my $newpos;

    if($nshared>0)
    {
        $newref=join("",splice(@aref,$nshared));
        $newalt=join("",splice(@aalt,$nshared));
        $newpos=$pos+$nshared;
    }
    else
    {
        $newref=$ref;
        $newalt=$alt;
        $newpos=$pos;
    }

    if($newref eq "")
    {
        $newref="-";
    }
    
    if($newalt eq "")
    {
        $newalt="-";
    }
    
    return(join($OFS,$chr,$newpos,$newref,$newalt));
}
