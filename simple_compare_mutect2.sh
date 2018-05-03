#!/bin/bash
module load bcftools/1.4.0
module load gatk/4.0.1.2

usage="Usage: $0 A.vcf B.vcf outname"

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]]
then
    echo -e $usage
    exit 1
fi

A=$1
B=$2
outname=$3

echo "SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity">$outname.csv
echo "SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity">${outname}_PASS.csv

#tail -n +2 pairs.csv | while read Patient A B Code Afile Bfile Nfile
#do

Afiltname=$(echo $A | sed "s/\.vcf/_filtered.vcf.gz/g")
Bfiltname=$(echo $B | sed "s/\.vcf/_filtered.vcf.gz/g")

gatk FilterMutectCalls -V $A -O $Afiltname
#--contamination-table tumor_calculatecontamination.table \
gatk FilterMutectCalls -V $B -O $Bfiltname

bcftools isec $Afiltname $Bfiltname -p ${outname}_PASS -f .,PASS
privA=$(sed -n '/^#.*/!p' ${outname}_PASS/0000.vcf | wc -l)
privB=$(sed -n '/^#.*/!p' ${outname}_PASS/0001.vcf | wc -l)
common=$(sed -n '/^#.*/!p' ${outname}_PASS/0002.vcf | wc -l)
sim=$(perl -e "print($common/($privA+$privB+$common))")
echo ${A},${B},$privA,$privB,$common,$sim >> ${outname}_PASS.csv

bcftools isec $Afiltname $Bfiltname -p ${outname}
privA=$(sed -n '/^#.*/!p' ${outname}/0000.vcf | wc -l)
privB=$(sed -n '/^#.*/!p' ${outname}/0001.vcf | wc -l)
common=$(sed -n '/^#.*/!p' ${outname}/0002.vcf | wc -l)
sim=$(perl -e "print($common/($privA+$privB+$common))")
echo ${A},${B},$privA,$privB,$common,$sim >> ${outname}.csv

#done
