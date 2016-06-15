#!/bin/bash
#

module load annovar/2014-07-14

humandb_dir=~/ngcchome/humandb/

outputfile=$(basename $1)
if [[ -f ${outputfile}.annotated.log ]]
then
    echo "This file had already been anotated"
else
    convert2annovar.pl --format vcf4 $1 > ${outputfile}.inputann
    annotate_variation.pl --geneanno --buildver hg19 --outfile ${outputfile}.annotated ${outputfile}.inputann $humandb_dir
    rm -f ${outputfile}.inputann
fi
