#!/bin/bash
#SBATCH -t 4-00:00

usage="$0 directory torun_file exe_params filtering_params NAB_params NAB2_params n_cores"

if [[ $# -ne 7 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]  || [[ ! -f $3 ]]  || [[ ! -f $4 ]]  || [[ ! -f $5 ]]  || [[ ! -f $6 ]]
then
    echo $usage
    exit
else
    dir=$(readlink -f $1)
    torun=$(readlink -f $2)
    exe_params=$(readlink -f $3)
    filtering_params=$(readlink -f $4)
    NAB_params=$(readlink -f $5)
    NAB2_params=$(readlink -f $6)
    n_cores=$7
fi

while read -r output normal a b
do
(time $SCRIPTSVCF_DIR/vcfFilteringTableV2_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --NABfilt_cond_inputfile2 $NAB2_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores $n_cores > ${output}.out ) &
done < $torun

wait
