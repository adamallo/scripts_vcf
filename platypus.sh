#!/bin/bash
#
#SBATCH -n 1

module load parallel/20140822 platypus/0.8.1 python/2.7.9

bamfiles=$1
output=$2
logfilename=$3
filterINDELS=$4

if [[ $# -ge 4 ]]
then
    shift 4
    time python /packages/6x/platypus/0.8.1/platypus callVariants --nCPU=$SLURM_JOB_CPUS_PER_NODE --bamFiles=$bamfiles --refFile=$HUMAN_GENOME --output=$output --logFileName=$logfilename $@ 
else
    echo "Error, the number of arguments for this script is not appropriate"
fi

mv $output ${output}_bkp_multisnv
$SCRIPTSVCF_DIR/perl.sh $SCRIPTSVCF_DIR/separateMultipleSNVPlatypus.pl -i ${output}_bkp_multisnv -o $output -f $filterINDELS
