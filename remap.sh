#!/bin/bash

module load tophat/2.0.14

usage="$0 refgenome refgenomenoy bamfile"

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ ! -f $3 ]]
then
    echo $usage
    exit
fi

name=$(basename $3 | sed "s/\.bam//")
dir=$(readlink -f $3 | sed "s/\/[^\/]*$/\//")
infile=$3
refgenome=$1
refgenomenoy=$2

#Reconstructs fastq from BAM
samtools sort -n $infile > $dir/${name}_sortedname.bam
bam2fastx --fastq --all -o $dir/${name}_reads.fq -N -P $dir/${name}_sortedname.bam
rm -rf $dir/${name}_sortedname.bam
#bamtools convert -in $infile -format fastq > $dir/${name}_reads.fq #Reconstructs fastq from BAM

bwa mem $refgenomenoy $dir/${name}_reads.1.fq $dir/${name}_reads.2.fq | samtools view -b > $dir/${name}_noY.bam #Remaps to noY
bwa mem $refgenome $dir/${name}_reads.1.fq $dir/${name}_reads.2.fq | samtools view -b > $dir/${name}_withY.bam #Remaps to normal
