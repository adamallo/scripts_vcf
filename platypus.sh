#!/bin/bash
#

cd $PBS_O_WORKDIR
module load parallel platypus
if [[ $# -eq 3 ]]
then
        time python /cm/shared/apps/platypus/Platypus.py callVariants --nCPU=$PBS_NUM_PPN --bamFiles=$1 --refFile=/home/afortun2/DCIS/GRCh37-lite.fa --output=$2 --logFileName=$3
elif [[ $# -ge 4 ]]
then
    bamfiles=$1
    output=$2
    logfilename=$3
    shift 3
    time python /cm/shared/apps/platypus/Platypus.py callVariants --nCPU=$PBS_NUM_PPN --bamFiles=$bamfiles --refFile=/home/afortun2/DCIS/GRCh37-lite.fa --output=$output --logFileName=$logfilename $@ 
else
    echo "Error, the number of arguments for this script is not appropriate"
fi
