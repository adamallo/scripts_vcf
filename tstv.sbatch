#!/bin/bash
#
#SBATCH -p private
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 4-00:00

module load vcftools/0.1.12b
module load parallel/20140822

if [[ $# -ne 1 ]] || [[ ! -d $1 ]]
then
    echo "$0 input_folder"
    exit
fi

dir=$1

cd $dir

files=$(ls -U *.vcf)

fun() {
    
    name=$(echo $1 | sed "s/.vcf//")
    ratio=$(vcftools --vcf $name.vcf --TsTv-summary --out $name 2>&1 | sed -n -e "/Ts\/Tv ratio/p" | sed "s/Ts\/Tv ratio: //")
    echo "$name,$ratio"
}

export -f fun

echo "File,tstv" > tstv.csv

parallel -j $SLURM_JOB_CPUS_PER_NODE --joblog tstv.log fun {1} ::: $files >> tstv.csv

#for vcf in $dir/*.vcf
#do
#    name=$(echo $vcf | sed "s/.vcf//")
#    ratio=$(vcftools --vcf $name.vcf --TsTv-summary --out $name 2>&1 | sed -n -e "/Ts\/Tv ratio/p" | sed "s/Ts\/Tv ratio: //")
#    echo "$name $ratio" >> $outputfile
#done
