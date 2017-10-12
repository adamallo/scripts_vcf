#!/bin/bash

usage="Usage: $0 seqzfileA seqzfileB outdir patient"

if [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ $# -ne 4 ]]
then
	echo $usage
	exit 1
fi


Rscript $SCRIPTSVCF_DIR/cnvcallandcomp.R $1 $2 $3 $SLURM_CPUS_PER_TASK $4
