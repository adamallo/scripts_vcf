#!/bin/bash
#

module load annovar/2014-07-14

if [ -z ${HUMANDB_DIR+x} ]
then
    echo "You need to use the environment variable HUMANDB_DIR to point to the human genome for annovar\n";
    exit
fi

outputfile=$(basename $1)
if [[ -f ${outputfile}.annotated.log ]]
then
    echo "This file had already been anotated"
else
    convert2annovar.pl --format vcf4 $1 > ${outputfile}.inputann
    annotate_variation.pl --geneanno --buildver hg19 --outfile ${outputfile}.annotated ${outputfile}.inputann $HUMANDB_DIR
    rm -f ${outputfile}.inputann
fi
