#!/bin/bash

usage="\n$0 directory torun_file exe_params filtering_params NAB_params NAB2_params covB_params popAF_params n_cores output_vcf output_list comprehensiveness [queue] \n\ntorun_file structure: output N_file A_file B_file\n-------------------------------------------------\n
\n
This script executes HeterAnalyzer_control.pl for each sample in a directory with its name. Then it integrates all the information in a file named results.csv and results_basictstv.csv\n"

queue=""

if (! ([[ $# -eq 12 ]] || [[ $# -eq 13 ]])) || [[ ! -d $1 ]] || [[ ! -f $2 ]]  || [[ ! -f $3 ]]  || [[ ! -f $4 ]]  || [[ ! -f $5 ]]  || [[ ! -f $6 ]] || [[ ! -f $7 ]] || [[ ! -f $8 ]]

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
    popAF_params=$(readlink -f $8)
    n_cores=$9
    output_vcf=${10}
    output_list=${11}
    comp=${12}
    queue=${13}
fi

EXE_DIR=$SCRIPTSVCF_DIR
#declare -A job_ids

dependency=""
while read -r output normal a b dnac
do
    if [[ $queue == "" ]]
    then
        id=$($EXE_DIR/HeterAnalyzer_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --NABfilt_cond_inputfile2 $NAB2_params --covaltB_cond_inputfile $covB_params --popAF_cond_inputfile $popAF_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores $n_cores --output_vcf $output_vcf --output_list $output_list --comp $comp | tee $dir/${output}.out | tail -n 1)
    else
        id=$($EXE_DIR/HeterAnalyzer_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --NABfilt_cond_inputfile2 $NAB2_params --covaltB_cond_inputfile $covB_params --popAF_cond_inputfile $popAF_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores $n_cores --output_vcf $output_vcf --output_list $output_list --comp $comp --queue $queue | tee $dir/${output}.out | tail -n 1) 
    fi
    dependency="${dependency}:${id}"
done < $torun

tstv=1

if [[ $comp -ne 2 ]] || [[ ${10} == 0 ]]
then
    tstv=0
fi

##Commented out
#if [[ $queue == "" ]]
#then
#    submit --dependency=afterok$dependency $EXE_DIR/postHeterAnalyzer.sh $dir $torun $exe_params $filtering_params $NAB_params $NAB2_params $covB_params $popAF_params $n_cores $tstv
#else
#    submit --partition=$queue --dependency=afterok$dependency $EXE_DIR/postHeterAnalyzer.sh $dir $torun $exe_params $filtering_params $NAB_params $NAB2_params $covB_params $popAF_params $n_cores $tstv
#fi 