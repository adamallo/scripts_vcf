##Output post-processing
for i in {B3,B6,D5,D8,DCIS64,K12}; do perl ~/ngcchome/dcis/scripts_vcf/tstv_pretabulate.pl ../$i/tstv.csv ../$i/tstv.tomerge; perl ~/ngcchome/dcis/scripts_vcf/merge_csv.pl ${i}.csv ../$i/tstv.tomerge ../${i}_withtstv.csv;done
for i in *.csv; do name=$(basename $i | sed "s/\.csv//"); ~/ngcchome/dcis/scripts_vcf/tabulate_results.pl -i $i -o ${name}_tab.csv --remove_all_but filtN_N,filtN_prop_mean,filtNAB_N,filtNAB_prop_mean,TsTv_A,TsTv_B,TsTv_filt,TsTv_filtN,TsTv_filtNAB ;done
for i in *; do cat $i | sed "s/[0-9]_/0_/g" > $i.mod;done
for i in *.csv; do mv $i $i.bkp;done
for i in *.mod; do name=$(echo $i | sed "s/.mod//"); mv $i $name;done

#Indexation of bam files
for i in *.bam; do sbatch -p private <(echo -e '#!/bin/bash'"\nmodule load samtools/1.3.1\n samtools index $i");done
