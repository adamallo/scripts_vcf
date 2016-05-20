#!/bin/bash
#

cd $PBS_O_WORKDIR
module load annovar

humandb_dir=~/humandb/

if [[ ! -d $1 ]]
then
	echo "Usage: $0 dir_with_vcfs"
	exit
fi

cd $1

for i in *filtN*.vcf
do 
    outputfile=$(basename $i)
    sem -j $PBS_NUM_PPN convert2annovar.pl --format vcf4 $1 > ${outputfile}.inputann ";" annotate_variation.pl --geneanno --buildver hg19 --outfile "${outputfile}.annotated" "${outputfile}.inputann" $humandb_dir
done
sem --wait
