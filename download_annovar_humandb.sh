#!/bin/bash
#

cd $PBS_O_WORKDIR
module load annovar

humandb_dir=~/humandb
annotate_variation.pl --downdb refGene $humandb_dir --build hg19
