#!/bin/bash

usage="$0 directory torun_file"

if [[ $# -ne 2 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]] 
then
    echo $usage
    exit
else
    dir=$1
    torun=$2
fi

while read -r output normal a b
do
(time ./vcfFilteringTableV2_control.pl -e $dir/exe_params -f $dir/filtering_params --NABfilt_cond_inputfile $dir/NAB_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores 16 > ${output}.out ) &
done < $torun

wait

