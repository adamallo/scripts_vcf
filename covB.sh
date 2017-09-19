#!/bin/bash

usage="$0 vcf1 vcf2 outputvcf_datavcf1bedvcf2 outputvcf_datavcf2bedvcf1 \n Environment variable GATKJAR should point to GATK\'s jar executable"

if [[ $# -ne 4 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit 1
fi

vcf1=$1
vcf2=$2
out1=$3
out2=$4

module load module load gatk/3.5.0

DATA=/home/dmalload/dcis/problem_data/new_data
GENOME=/home/dmalload/my_storage/GRCh37-lite.f

dothething ()
{
    vcf1=$1
    vcf2=$2
    out=$3

    if [[ ! -f $out ]]
    then
        name_vcf2=$(echo $2 | sed "s/.vcf//g")
        name_out=$(echo $out | sed "s/.tsv//g")
        vcf2bed --deletions < $vcf2 > ${name_vcf2}_deletions.bed
        vcf2bed --snvs < $vcf2 > ${name_vcf2}_snvs.bed
        bedops --everything ${name_vcf2}_{deletions,snvs}.bed | awk 'BEGIN{OFS="\t"}{print($1,$2,$3)}' > ${name_vcf2}.bed
        ###The input for GATK is wrong
        java -Xms512m -Xmx6G -jar $GATKJAR -T UnifiedGenotyper -R $GENOME -I $vcf1 -o "$name_out.vcf" --intervals ${name_vcf2}.bed --output_mode EMIT_ALL_SITES > "$name_out.log" 2>&1
        cat "$name_out.vcf" | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[2]=="" and $reads[1]=0;if($reads[0] eq "./."){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(",",@reads)};print join("\t",@F[0,1,3,4],$readsref,$readsout)' > $out
    fi
}

dothething $vcf1 $vcf2 $out1
dothething $vcf2 $vcf1 $out2
