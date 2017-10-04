#!/bin/bash

usage="Usage: $0 seqzfile outdir"

if [[ ! -f $1 ]] || [[ $# -ne 2 ]]
then
	echo $usage
	exit 1
fi


Rscript $BLADDER_SCRIPTS/sequenza.R $1 $2 $SLURM_CPUS_PER_TASK
