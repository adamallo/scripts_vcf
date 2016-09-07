#!/bin/bash
#
#SBATCH -n 1

module load parallel/20140822 platypus/0.8.1 python/2.7.9

if [[ $# -eq 3 ]]
then
        time python /packages/6x/platypus/0.8.1/platypus callVariants --nCPU=$SLURM_JOB_CPUS_PER_NODE --bamFiles=$1 --refFile=/home/dmalload/temp_storage/GRCh37-lite.fa --output=$2 --logFileName=$3
elif [[ $# -ge 4 ]]
then
    bamfiles=$1
    output=$2
    logfilename=$3
    shift 3
    time python /packages/6x/platypus/0.8.1/platypus callVariants --nCPU=$SLURM_JOB_CPUS_PER_NODE --bamFiles=$bamfiles --refFile=/home/dmalload/temp_storage/GRCh37-lite.fa --output=$output --logFileName=$logfilename $@ 
else
    echo "Error, the number of arguments for this script is not appropriate"
fi
