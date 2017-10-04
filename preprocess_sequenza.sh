#!/bin/bash
#SBATCH --mem 10000

module unload python
module load pypy27/5.6

usage="Usage: $0 tumor normal outdir outname"

if [[ $# -ne 4 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit 1
fi

gcfile="$(dirname $HUMAN_GENOME)/b37.gc50Base.txt.gz"

pypy `which sequenza-utils.py` pileup2seqz -gc $gcfile -t $1 -n $2 | gzip > $3/$4_temp.seqz.gz
#pypy `which sequenza-utils.py` bam2seqz -gc $gcfile --fasta $HUMAN_GENOME -t $1 -n $2 | gzip > $3/$4_temp.seqz.gz
pypy `which sequenza-utils.py` seqz-binning -w 50 -s $3/$4_temp.seqz.gz | gzip > $3/$4.seqz.gz
#rm -f $3/$4_temp.seqz.gz
