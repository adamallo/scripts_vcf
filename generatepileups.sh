#!/bin/bash
#SBATCH --mem 10000

usage="Usage: $0 outputfile bamfile"

if [[ ! -f $2 ]] || [[ $# -ne 2 ]]
then
    echo -e $usage
    exit 1
fi

##Filtering alignments with mapQ<10 and bases with Qual<20
samtools mpileup -f $HUMAN_GENOME -q 10 -Q 20 $2 | gzip > $1
