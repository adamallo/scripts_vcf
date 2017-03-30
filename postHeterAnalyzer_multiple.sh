#!/bin/bash

usage="\n$0 directory torun_file exe_params filtering_params NAB_params NAB2_params\n\ntorun_file structure: output N_file A_file B_file\n-------------------------------------------------\n
\n
This script postprocess the output of HeterAnalyzer_control.pl for each sample in a directory with its name, integrating all the information in a file named results.csv and results_basictstv.csv\n"

if [[ $# -ne 6 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]  || [[ ! -f $3 ]]  || [[ ! -f $4 ]]  || [[ ! -f $5 ]]  || [[ ! -f $6 ]]
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
fi

flag=0

EXE_DIR=$SCRIPTSVCF_DIR

while read -r output normal a b
do
    ##General tstv postprocessing
    i=$dir/$output
    namenormal=$(basename $normal )
    namenormal=$(echo $namenormal | sed "s/.bam//")
    namea=$(basename $a)
    namea=$(echo $namea | sed "s/.bam//")
    nameb=$(basename $b)
    nameb=$(echo $nameb | sed "s/.bam//")

    tail -n +2 tstv.csv | sed -n -e "/^${namenormal},/p" -e "/^${namea}.*,/p" -e "/^${nameb}.*,/p" -e "/^${namenormal}_${namea}_${nameb}.*,/p" -e "/^${namenormal}_${nameb}_${namea}.*,/p" | sed -e "s/$namenormal/N/g" -e "s/$namea/A/g" -e "s/$nameb/B/g" -e "s/N_A_B/NAB/" -e "s/N_B_A/NAB/" | sort | uniq >> $i/tstv.csv
    perl $EXE_DIR/tstv_pretabulate_heter.pl $i/tstv.csv $i/tstv.tomerge
    #name=$(basename $i)
    #dir=$(echo $i | sed "s/\/[^\/]*$//")
    cat ${i}/${output}.csv | awk 'BEGIN{OFS=",";FS=","}{out="";for (i=2;i<=NF;++i){out=out $i OFS};out=substr(out,0,length(out)-1);print out}' > $i/csvtomerge.temp
    perl $EXE_DIR/merge_csv.pl $i/csvtomerge.temp $i/tstv.tomerge $i/${output}_withtstv.temp
    rm -f $i/csvtomerge.temp
    cat $i/${output}_withtstv.temp | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Sample\",\$0} else {print \"$output\",\$0}}" > ${dir}/${output}_withtstv.csv
    rm -f $i/${output}_withtstv.temp

    #BasicTsTv postprocessing and general tstv data gathering
    a=$(cat $i/tstv.csv | sed -n "/^A,/p" | sed "s/^A,//")
    b=$(cat $i/tstv.csv | sed -n "/^B,/p" | sed "s/^B,//")
    n=$(cat $i/tstv.csv | sed -n "/^N,/p" | sed "s/^N,//")

    if [ $flag -ne 0 ]
    then
        tail -n +2 $dir/${output}_withtstv.csv >> $dir/results.csv
    else
        cat $dir/${output}_withtstv.csv > $dir/results.csv
        echo "Sample,TsTv_A,TsTv_B,TsTv_N" > $dir/results_basictstv.csv
        flag=1
    fi

    echo "$output,$a,$b,$n" >> $dir/results_basictstv.csv

done < $torun
