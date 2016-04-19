#!/bin/bash
#
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR
module load parallel platypus
if [[ $# -eq 3 ]]
then
        time python /cm/shared/apps/platypus/Platypus.py callVariants --nCPU=8 --bamFiles=$1 --refFile=/home/afortun2/DCIS/GRCh37-lite.fa --output=$2 --logFileName=$3

elif [[ $# -eq 4 ]]
then
        time python /cm/shared/apps/platypus/Platypus.py callVariants --nCPU=8 --bamFiles=$1 --refFile=/home/afortun2/DCIS/GRCh37-lite.fa $2 --output=$3 --logFileName=$4
else
    echo "Error, the number of arguments for this script is not appropriate"
fi
