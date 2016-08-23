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
(time $SCRIPTSVCF_DIR/vcfFilteringTableV2_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --NABfilt_cond_inputfile2 $NAB2_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores $n_cores > $dir/${output}.out ) &
done < $torun

wait

flag=0

while read -r output normal a b
do
    i=$dir/$output
    perl $SCRIPTSVCF_DIR/tstv_pretabulate.pl $i/tstv.csv $i/tstv.tomerge
    perl $SCRIPTSVCF_DIR/merge_csv.pl ${i}.csv $i/tstv.tomerge ${i}_withtstv.csv
    sed -i.bak "s/[0-9]_/0_/g" ${i}_withtstv.csv ## The normal id may change across replicates but it is always the same condition. We change them all to 0.
    perl $SCRIPTSVCF_DIR/tabulate_results.pl -i ${i}_withtstv.csv -o ${i}_tab.csv 
    cat ${i}_tab.csv | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Sample\",\$0} else {print \"$output\",\$0}}" > ${dir}/${output}_tab_sample.csv

    #BasicTsTv postprocessing and general tstv data gathering
    a=$(cat $i/tstv.csv | sed -n "/^A,/p" | sed "s/^A,//")
    b=$(cat $i/tstv.csv | sed -n "/^B,/p" | sed "s/^B,//")
    n=$(cat $i/tstv.csv | sed -n "/^N,/p" | sed "s/^N,//")

    if [ $flag -ne 0 ]
    then
        tail -n +2 $dir/${output}_tab_sample.csv >> $dir/results.csv
    else
        cat $dir/${output}_tab_sample.csv > $dir/results.csv
        echo "Sample,TsTv_A,TsTv_B,TsTv_N" > $dir/results_basictstv.csv
        flag=1
    fi

    echo "$output,$a,$b,$n" >> $dir/results_basictstv.csv
done < $torun
