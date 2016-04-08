#!/bin/bash
#
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR
module load parallel platypus
time python /cm/shared/apps/platypus/Platypus.py callVariants --nCPU=8 --bamFiles=$1 --refFile=/home/afortun2/DCIS/GRCh37-lite.fa $2 --output=$3 --logFileName=$4
