#!/bin/bash

usage="$0 outputname vcf1 vcf2 ... vcfn\n Environment variable GATKJAR should point to GATK\'s jar executable"

if [[ $# -lt 2 ]]
then
    echo -e $usage
    exit 1
fi

module load bedops/2.4.35

out=$1

shift

DATA=/home/dmalload/dcis/problem_data/new_data
GENOME=/home/dmalload/my_storage/GRCh37-lite.fa

beds=""

if [[ ! -f $out ]]
then
    for vcf in "$@"
    do
        if [[ ! -s $vcf ]]
        then
            echo "Error reading the file $vcf"
            echo -e $usage
            exit 1
        fi
        vcf2bed --deletions < $vcf > ${vcf}_deletions.bed
        vcf2bed --insertions < $vcf > ${vcf}_insertions.bed
        vcf2bed --snvs < $vcf > ${vcf}_snvs.bed
        beds="$beds${vcf}_deletions.bed ${vcf}_snvs.bed ${vcf}_insertions.bed "
    done
    ##working here paste $@ comma separated into allvcfs
    bedops --everything $beds | awk 'BEGIN{OFS="\t"}{print($1,$2,$3)}' > $out
else
    echo "The file $out is present and will be reused"
fi
