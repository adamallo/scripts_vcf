#!/bin/bash
#SBATCH --mem=10G

usage="$0 input1 input output metrics\n"

if [[ $# -ne 4 ]] || [[ ! -s $1 ]] || [[ ! -s $2 ]]
then
    echo -e $usage
    exit 1
fi

module load picard/2.9.2
bam1=$1
bam2=$2
finalbam=$3
metrics=$4
tempbam=$(echo $finalbam | sed "s/.bam/_temp.bam/")
tempdir="/scratch/dmalload/tmp/$(echo $finalbam | sed "s/\.[b|s]am//I")"
mkdir -p $tempdir

if [[ ! -s $tempbam ]]
then
    java -Xmx8G -jar $PICARDJAR MergeSamFiles I=$bam1 I=$bam2 O=$tempbam TMP_DIR=$tempdir SORT_ORDER=queryname
fi

java -Xmx8G -jar $PICARDJAR MarkDuplicates TMP_DIR=$tempdir I=$tempbam O=/dev/stdout M=$metrics CREATE_INDEX=false | sambamba sort --tmpdir=$tempdir -o $finalbam /dev/stdin

if [[ -s $finalbam ]] && [[ $? -eq 0 ]]
then
    rm -rf $tempbam
fi

rm -rf $tempdir
