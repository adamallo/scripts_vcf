#!/bin/bash
#SBATCH --mem 10000

##Loading my sequenza environment
module unload python
module load pypy27/5.6
source ~/virtualenvs/sequenza/bin/activate

usage="Usage: $0 tumorBam normalBam outdir outname\nIf pileup.gz files are available, this script will use them instead of the bam files for an speed increase."

if [[ $# -ne 4 ]]
then
    echo -e $usage
    exit 1
fi

tumor=$1
normal=$2
outdir=$3
outname=$4

gcfile="$(dirname $HUMAN_GENOME)/b37.gc50Base.txt.gz"

if [[ ! -f $gcfile ]]
then
    sequenza-utils gc_wiggle --fasta $HUMAN_GENOME -w 50 -o $gcfile
fi

pileupTumor=$(echo $tumor | sed "s/.bam/.pileup.gz/")
pileupNormal=$(echo $normal | sed "s/.bam/.pileup.gz/")

if [[ -f $pileupTumor ]] && [[ -f $pileupNormal ]] 
then
    sequenza-utils bam2seqz -p -gc $gcfile -n $pileupNormal -t $pileupTumor --fasta $HUMAN_GENOME -o $outdir/${outname}_temp.seqz.gz
elif [[ -f $tumor ]] && [[ -f $normal ]]
then
    ##Default, using bam files. It will involve the conversion of bam files to pileups, so it is more efficient to first do that independently for each file. Especially for multi-sample experiments.
    sequenza-utils bam2seqz -gc $gcfile --fasta $HUMAN_GENOME -t $tumor -n $normal -o $outdir/${outname}_temp.seqz.gz
else
    echo ERROR: $normal and/or $tumor not found
fi 

#Size reduction
sequenza-utils seqz_binning -w 50 -s $outdir/${outname}_temp.seqz.gz -o $outdir/${outname}.seqz.gz

if [[ $? -eq 0 ]]
then
    rm -f $outdir/${outname}_temp.seqz.gz $outdir/${outname}_temp.seqz.gz.tbi
fi
