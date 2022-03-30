#!/bin/bash

usage="Usage: $0 seqzfile outdir"

if [[ ! -f $1 ]] || [[ $# -ne 2 ]]
then
	echo $usage
	exit 1
fi

module unload gcc/4.9.2
module load r/4.1.2-BLAS

Rscript $SCRIPTSVCF_DIR/sequenza.R $1 $2 $SLURM_CPUS_PER_TASK
