#!/bin/bash
#

cd $PBS_O_WORKDIR
module load annovar

humandb_dir=~/humandb/

outputfile=$(basename $1)
convert2annovar.pl --format vcf4 $1 > ${outputfile}.inputann
annotate_variation.pl --geneanno --buildver hg19 --outfile "${outputfile}.annotated" "${outputfile}.inputann" $humandb_dir
