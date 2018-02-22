#!/bin/bash
#SBATCH -c 8

#usage="Usage: $0 bamfile"
usage="Usage: $0 commonname"

#if [[ ! -f $1 ]]
if [[ ! -f ${1}_R1.fastq.gz ]] || [[ ! -f ${1}_R2.fastq.gz ]]
then
    echo -e $usage
    exit 1
fi

#skewer -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
#    -Q 10 -z -t $SLURM_CPUS_PER_TASK -n ${1}_R1.fastq.gz ${1}_R2.fastq.gz -o ${1}
skewer -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -Q 10 -z -t $SLURM_CPUS_PER_TASK -n ${1}_R1.fastq.gz ${1}_R2.fastq.gz -o ${1}
