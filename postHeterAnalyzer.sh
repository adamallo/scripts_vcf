#!/bin/bash

usage="\n$0 directory torun_file exe_params filtering_params NAB_params NAB2_params covB_params n_cores tstv\n\ntorun_file structure: output N_file A_file B_file\n-------------------------------------------------\n
\n
This script postprocess the output of HeterAnalyzer_control.pl for each sample in a directory with its name, integrating all the information in a file named results.csv. Results results_basictstv.csv is generated only if tstv==1.\n"

if [[ $# -ne 9 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]  || [[ ! -f $3 ]]  || [[ ! -f $4 ]]  || [[ ! -f $5 ]]  || [[ ! -f $6 ]] || [[ ! -f $7 ]]
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
    tstv=$9
fi

flag=0

EXE_DIR=$SCRIPTSVCF_DIR

while read -r output normal a b
do
    ##General tstv postprocessing
    i=$dir/$output
    #name=$(basename $i)
    #dir=$(echo $i | sed "s/\/[^\/]*$//")
    cat ${dir}/${output}.csv | awk 'BEGIN{OFS=",";FS=","}{out="";for (i=2;i<=NF;++i){out=out $i OFS};out=substr(out,0,length(out)-1);print out}' > $i/csvtomerge.temp
    
    if [[ $tstv -eq 1 ]]
    then
        perl $EXE_DIR/tstv_pretabulate_heter.pl $i/tstv.csv $i/tstv.tomerge
        perl $EXE_DIR/merge_csv.pl $i/csvtomerge.temp $i/tstv.tomerge $i/${output}_withtstv.temp
        rm -f $i/csvtomerge.temp
        perl $SCRIPTSVCF_DIR/tabulate_results.pl -i $i/${output}_withtstv.temp -o $i/${output}_withtstv_tab.temp
        cat $i/${output}_withtstv_tab.temp | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Sample\",\$0} else {print \"$output\",\$0}}" > ${dir}/${output}_withtstv.csv
        rm -f $i/${output}_withtstv.temp $i/${output}_withtstv_tab.temp

        #BasicTsTv postprocessing and general tstv data gathering
        a=$(cat $i/tstv.csv | sed -n "/^A,/p" | sed "s/^A,//")
        b=$(cat $i/tstv.csv | sed -n "/^B,/p" | sed "s/^B,//")
        n=$(cat $i/tstv.csv | sed -n "/^N,/p" | sed "s/^N,//")
    else 
        perl $SCRIPTSVCF_DIR/tabulate_results.pl -i $i/csvtomerge.temp -o $i/csvtomerge_tab.temp
        cat $i/csvtomerge_tab.temp | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Sample\",\$0} else {print \"$output\",\$0}}" > ${dir}/${output}_withtstv.csv
        rm -f $i/csvtomerge.temp $i/csvtomerge_tab.temp
    fi
    if [ $flag -ne 0 ]
    then
        tail -n +2 $dir/${output}_withtstv.csv >> $dir/results.csv
    else
        cat $dir/${output}_withtstv.csv > $dir/results.csv
        if [[ $tstv -eq 1 ]]
        then
            echo "Sample,TsTv_A,TsTv_B,TsTv_N" > $dir/results_basictstv.csv
        fi
        flag=1
    fi
    if [[ $tstv -eq 1 ]]
    then
        echo "$output,$a,$b,$n" >> $dir/results_basictstv.csv
    else
        rm -f ${dir}/${output}_withtstv.csv
    fi
done < $torun
