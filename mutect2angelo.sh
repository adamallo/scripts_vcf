#!/bin/bash
#SBATCH --mem 10000
module load samtools/1.4.0
module load gatk/4.0.1.2

usage="Usage: $0 normal tumor outname"

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit 1
fi

output=$3
tumor=$2
normal=$1

tumorname=$(samtools view -H $tumor | sed -n '/@RG/p' | sed "s/.*SM:\([^ \t]*\).*/\1/g");
if [[ $tumorname -eq "" ]]
then
    echo "Error, the tumor bam file $tumor does not have a readgroup"
    exit
fi
gatk Mutect2 -R $HUMAN_GENOME -I $normal -I $tumor --tumor-sample $tumorname --output $output.vcf
