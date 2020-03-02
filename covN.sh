#!/bin/bash

usage="$0 vcf1_name vcf2_name outputname Nbam\n Environment variable GATKJAR should point to GATK\'s jar executable"

if [[ $# -ne 4 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ ! -f $4 ]]
then
    echo -e $usage
    exit 1
fi

vcf1=$1
vcf2=$2
out=$3
bam1=$4

module load gatk/3.5.0
module load bedops/2.4.35

DATA=/home/dmalload/dcis/problem_data/new_data
GENOME=/home/dmalload/my_storage/GRCh37-lite.fa

if [[ ! -f $out ]]
then
    name_vcf1=$(echo $vcf1 | sed "s/.vcf//g")
    name_vcf2=$(echo $vcf2 | sed "s/.vcf//g")
    name_out=$(echo $out | sed "s/.tsv//g")
    vcf2bed --deletions < $vcf2 > ${name_vcf2}_covN_deletions.bed
    vcf2bed --insertions < $vcf2 > ${name_vcf2}_covN_insertions.bed
    vcf2bed --snvs < $vcf2 > ${name_vcf2}_covN_snvs.bed
    vcf2bed --deletions < $vcf1 > ${name_vcf1}_covN_deletions.bed
    vcf2bed --insertions < $vcf1 > ${name_vcf1}_covN_insertions.bed
    vcf2bed --snvs < $vcf1 > ${name_vcf1}_covN_snvs.bed
    bedops --everything {${name_vcf2},${name_vcf1}}_covN_{deletions,insertions,snvs}.bed | awk 'BEGIN{OFS="\t"}{print($1,$2,$3)}' > ${name_out}_covN.bed
    java -Xms512m -Xmx6G -jar $GATKJAR -T UnifiedGenotyper -R $GENOME -I $bam1 -o "$name_out.vcf" --intervals ${name_out}_covN.bed --output_mode EMIT_ALL_SITES -glm BOTH > "$name_out.log" 2>&1
    cat "$name_out.vcf" | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[1]=="" and $reads[1]=0;if($reads[0] eq "./."){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(",",@reads)};print join("\t",@F[0,1,3,4],$readsref,$readsout)' > $out
else
    echo "The file $out is present and will be reused"
fi
