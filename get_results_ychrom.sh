#!/bin/bash

usage="$0 directory torun_file min_mapq max_mapq output_file"

if [[ $# -ne 5 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]
then
    echo $usage
    exit
else
    dir=$(readlink -f $1)
    torun=$(readlink -f $2)
    min_mapq=$3
    max_mapq=$4
    output_file=$5
fi

echo "Patient,Sample,DCIS,Min_mapq,Max_mapq,Total,CHRY,CHRY-NOY,CHRY-NONPAR,CHRY-NONPAR-NOY,CHRY-NONPAR-BED,CHRY-NONPAR-BED-NOY" > $output_file
while read -r output normal a b
do
patient=$(echo $output | sed 's/_.\.out//')
awk -v patient="${patient}" -v dcis="0" -v sample="N" -v mapq="${min_mapq}" -v max_mapq="${max_mapq}" 'BEGIN{FS=",";OFS=","}{print patient,sample,dcis,mapq,max_mapq,$1,$2,$3,$4,$5,$6,$7}' $dir/${output}_N.out >> $output_file
awk -v patient="${patient}" -v dcis="1" -v sample="A" -v mapq="${min_mapq}" -v max_mapq="${max_mapq}" 'BEGIN{FS=",";OFS=","}{print patient,sample,dcis,mapq,max_mapq,$1,$2,$3,$4,$5,$6,$7}' $dir/${output}_A.out >> $output_file
awk -v patient="${patient}" -v dcis="1" -v sample="B" -v mapq="${min_mapq}" -v max_mapq="${max_mapq}" 'BEGIN{FS=",";OFS=","}{print patient,sample,dcis,mapq,max_mapq,$1,$2,$3,$4,$5,$6,$7}' $dir/${output}_B.out >> $output_file
done < $torun
