while read -r output normal a b
do
(time ./vcfFilteringTableV2_control.pl -e second_round_params/exe_params -f second_round_params/filtering_params --NABfilt_cond_inputfile second_round_params/NAB_params -o ${output}.csv --normal_bamfile $normal --sample_A_bamfile $a --sample_B_bamfile $b --output_dir $output --n_cores 16 > ${output}.out ) &
done < $1
