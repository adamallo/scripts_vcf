#!/bin/bash
#SBATCH -c 8

#usage="Usage: $0 bamfile"
usage="Usage: $0 commonname"

#if [[ ! -f $1 ]]

pair1=$1-pair1.fastq.gz
pair2=$1-pair2.fastq.gz
gname=$1

if [[ ! -f $pair1 ]] || [[ ! -f $pair2 ]]
then
    echo -e $usage
    exit 1
fi

name=$(echo $gname | sed "s/.\+\/\(.\+\)-trimmed/\\1/g")
dirname=$(dirname $gname)

##No readgroup so far
#bamname=$(echo $gname | sed "s/\(.\+\)\/.\+/\\1.bam/g")

#name=$(echo $1 | sed "s/.bam//g")
#readgroup=$(samtools view -H $bamname | sed -n "/@RG.*$name.*/p" | sed "s/\t/\\t/g")

#bwa mem -R "$readgroup" -t $SLURM_CPUS_PER_TASK -M $HUMAN_GENOME "${dirname}/${name}-trimmed-pair1.fastq.gz" "${dirname}/${name}-trimmed-pair2.fastq.gz" | samtools sort -@ $SLURM_CPUS_PER_TASK -o "$dirname/${name}.bam"

bwa mem -t $SLURM_CPUS_PER_TASK -M $HUMAN_GENOME "${dirname}/${name}-trimmed-pair1.fastq.gz" "${dirname}/${name}-trimmed-pair2.fastq.gz" | samtools sort -@ $SLURM_CPUS_PER_TASK -o "$dirname/${name}.bam"

#cutadatp
#bwa mem -t $SLURM_CPUS_PER_TASK -M $HUMAN_GENOME ${name}_1_trimmed.fastq.gz ${name}_2_trimmed.fastq.gz | samtools sort -@ $SLURM_CPUS_PER_TASK -o ${name}_trimmed.bam
samtools flagstat "$dirname/${name}.bam" > "$dirname/${name}.samtools.txt"
