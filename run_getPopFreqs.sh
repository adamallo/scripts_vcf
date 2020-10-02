#!/bin/bash


for i in $1/*/vcfdict.csv
do 
    read a b c <<< $(perl -lne 'if(/AfiltcovBNAB.*different.vcf/){s/.*,//;$a=$_}elsif(/BfiltcovBNAB.*different.vcf/){s/.*,//;$b=$_}elsif(/filtcovBNABU.*common.vcf/){s/.*,//;$c=$_};END{print("$a $b $c")}' $i);
    const=$(readlink -e $i | sed 's/[^\/]\+$//') 
    files+=(a,$const$a b,$const$b c,$const$c);
    
done

runthething()
{
    $SCRIPTSVCF_DIR/perl.sh $SCRIPTSVCF_DIR/getPopFreqs.pl -p /home/dmalload/my_storage/gnomAD/gnomad.genomes.r2.1.sites.vcf.bgz -i $(echo $1 | sed 's/^.*,//') -o $(echo $1| sed 's/^\(.\),\(.*\)\/[^\/]\+$/\2\/\1_withPAF.tsv/')
}

export -f runthething

parallel -j $SLURM_JOB_CPUS_PER_NODE runthething {} ::: ${files[@]}
