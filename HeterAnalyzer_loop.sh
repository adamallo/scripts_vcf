#!/bin/bash

usage="$0 directory torun_file\n
torun_file structure: output N_file A_file B_file\n
-------------------------------------------------\n
\n
This script executes HeterAnalyzer_control.pl for each sample in a directory with its name. Then it integrates all the information in a file named results.csv\n"

if [[ $# -ne 2 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]] 
then
    echo -e $usage
    exit
else
    dir=$1
    torun=$2
fi

EXE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while read -r output normal a b
do
(time $EXE_DIR/HeterAnalyzer_control.pl -e $dir/exe_params -f $dir/filtering_params --NABfilt_cond_inputfile $dir/NAB_params --NABfilt_cond_inputfile2 $dir/NAB_params2 -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores 16 > $dir/${output}.out ) &
done < $torun

wait

flag=0

while read -r output normal a b
do
    if [ $flag -ne 0 ]
    then
        tail -n1 $dir/${output}.csv >> $dir/results.csv
    else
        cat $dir/${output}.csv > $dir/results.csv
        flag=1
    fi

done < $torun

