#!/bin/bash

usage="\n$0 directory torun_file exe_params filtering_params NAB_params NAB2_params n_cores\n\ntorun_file structure: output N_file A_file B_file\n-------------------------------------------------\n
\n
This script executes HeterAnalyzer_control.pl for each sample in a directory with its name. Then it integrates all the information in a file named results.csv and results_basictstv.csv\n"

if [[ $# -ne 7 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]  || [[ ! -f $3 ]]  || [[ ! -f $4 ]]  || [[ ! -f $5 ]]  || [[ ! -f $6 ]]
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
    n_cores=$7
fi

EXE_DIR=$SCRIPTSVCF_DIR
#declare -A job_ids

dependency=""
while read -r output normal a b
do
    id=$(sbatch -p private --job-name=${output}_HetAn <( echo -e '#!/bin/bash'"\n $EXE_DIR/HeterAnalyzer_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --NABfilt_cond_inputfile2 $NAB2_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores $n_cores > $dir/${output}.out") | sed "s/Submitted batch job \(.*\)/\1/")
    #(time $EXE_DIR/HeterAnalyzer_control.pl -e $exe_params -f $filtering_params --NABfilt_cond_inputfile $NAB_params --NABfilt_cond_inputfile2 $NAB2_params -o $dir/${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $dir/$output --n_cores $n_cores > $dir/${output}.out) &
    dependency="${dependency}:${id}"
    #job_ids[$id]=1
done < $torun

#dependency=${dependency:0:${#dependency} - 1}

#while [[ "${#job_ids[@]}" -ne 0 ]]
#do
#    for id in "${!job_ids[@]}"
#    do
#        qstat=$(qstat $id | tail -n 1 | awk 'BEGIN{FS=" ";var=1}{if ($5 != "C"){var=0}}END{print var}')
#        if [[ $qstat == 1 ]]
#        then
#            unset $job_ids[$id]
#        fi
#    done
#    sleep 10
#done

sbatch -p private --job-name=getdependencies --dependency=afterok$dependency <(echo -e '#!/bin/bash' "\ndependency=\"\";while read -r output normal a b;do id=\$(tail -n 1 $dir/\${output}.out);dependency=\"\$dependency:\$id\";done < $torun;sbatch -p private --dependency=afterok\$dependency $EXE_DIR/postHeterAnalyzer.sh $dir $torun $exe_params $filtering_params $NAB_params $NAB2_params $n_cores")

#TsTv calculation has been integrated in HeterAnalyzer_control

#flag=0
#
#while read -r output normal a b
#do
#    ##General tstv postprocessing
#    i=$dir/$output
#    perl $EXE_DIR/tstv_pretabulate_heter.pl $i/tstv.csv $i/tstv.tomerge
#    #name=$(basename $i)
#    #dir=$(echo $i | sed "s/\/[^\/]*$//")
#    cat ${dir}/${output}.csv | awk 'BEGIN{OFS=",";FS=","}{out="";for (i=2;i<=NF;++i){out=out $i OFS};out=substr(out,0,length(out)-1);print out}' > $i/csvtomerge.temp
#    perl $EXE_DIR/merge_csv.pl $i/csvtomerge.temp $i/tstv.tomerge $i/${output}_withtstv.temp
#    rm -f $i/csvtomerge.temp
#    cat $i/${output}_withtstv.temp | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Sample\",\$0} else {print \"$output\",\$0}}" > ${dir}/${output}_withtstv.csv
#    rm -f $i/${output}_withtstv.temp
#    
#    #BasicTsTv postprocessing and general tstv data gathering
#    a=$(cat $i/tstv.csv | sed -n "/^A,/p" | sed "s/^A,//")
#    b=$(cat $i/tstv.csv | sed -n "/^B,/p" | sed "s/^B,//")
#    n=$(cat $i/tstv.csv | sed -n "/^N,/p" | sed "s/^N,//")
#
#    if [ $flag -ne 0 ]
#    then
#        tail -n +2 $dir/${output}_withtstv.csv >> $dir/results.csv
#    else
#        cat $dir/${output}_withtstv.csv > $dir/results.csv
#        echo "Sample,TsTv_A,TsTv_B,TsTv_N" > $dir/results_basictstv.csv
#        flag=1
#    fi
#    
#    echo "$output,$a,$b,$n" >> $dir/results_basictstv.csv
#
#done < $torun

