#!/bin/bash
#SBATCH --mem 10000
module load picard/2.3.0

name=$(echo $1 | sed "s/.bam//g")

##Fix for a former problem using bwa mem. Now fixed and not necessary anymore

#samtools view -H ${name}.bam | sed -e "/@PG/s/\t/\\\t/g" -e "/@PG/s/\\\tID:\([^\\]\+\)\\\tPN:\([^\\]\+\)\\\tVN:\([^\\]\+\)\\\tCL:/\tID:\\1\tPN:\\2\tVN:\\3\tCL:/" | samtools reheader - ${name}.bam > ${name}_fixed.bam
#rm -f ${name}.bam
#mv ${name}_fixed.bam ${name}.bam

picard MarkDuplicates I=${name}.bam O=${name}_mdups.bam ASSUME_SORTED=true M=${name}.metrics
samtools index ${name}_mdups.bam
