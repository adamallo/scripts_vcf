#!/bin/bash

usage="\n$0 directory torun_file exe_params filtering_params NAB_params NAB2_params covB_params n_cores\n\ntorun_file structure: output N_file A_file B_file\n-------------------------------------------------\n
\n
This script executes HeterAnalyzer_control.pl for each sample in a directory with its name. Then it integrates all the information in a file named results.csv and results_basictstv.csv\n"

if [[ $# -ne 8 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]  || [[ ! -f $3 ]]  || [[ ! -f $4 ]]  || [[ ! -f $5 ]]  || [[ ! -f $6 ]] || [[ ! -f $7 ]]
then
    echo -e $usage
    exit
else
    dir=$(readlink -f $1)
    torun=$(readlink -f $2)
    exe_params=$(readlink -f $3)
    filtering_params=$(readlink -f $4)
    NAB_params=$(readlink -f $5)
    NAB2_params=$(readlink -f $6)
    covB_params=$(readlink -f $7)
    n_cores=$8
fi

EXE_DIR=$SCRIPTSVCF_DIR
#declare -A job_ids

dependency=""
while read -r output normal a b
do
    id=$($EXE_DIR/HeterAnalyzer_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --NABfilt_cond_inputfile2 $NAB2_params --covaltB_cond_inputfile $covB_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores $n_cores | tee $dir/${output}.out | tail -n 1)
    dependency="${dependency}:${id}"
done < $torun

submit --dependency=afterok$dependency $EXE_DIR/postHeterAnalyzer.sh $dir $torun $exe_params $filtering_params $NAB_params $NAB2_params $covB_params $n_cores

#sbatch -p private --job-name=getdependencies --dependency=afterok$dependency <(echo -e '#!/bin/bash' "\ndependency=\"\";while read -r output normal a b;do id=\$(tail -n 1 $dir/\${output}.out);dependency=\"\$dependency:\$id\";done < $torun;sbatch -p private --dependency=afterok\$dependency $EXE_DIR/postHeterAnalyzer.sh $dir $torun $exe_params $filtering_params $NAB_params $NAB2_params $covB_params $n_cores")
