##Output post-processing (replicated samples) (outdated
#for i in {B3,B6,D5,D8,DCIS64,K12}; do perl ~/ngcchome/dcis/scripts_vcf/tstv_pretabulate.pl ../$i/tstv.csv ../$i/tstv.tomerge; perl ~/ngcchome/dcis/scripts_vcf/merge_csv.pl ${i}.csv ../$i/tstv.tomerge ../${i}_withtstv.csv;done
#for i in *.csv; do name=$(basename $i | sed "s/\.csv//"); ~/ngcchome/dcis/scripts_vcf/tabulate_results.pl -i $i -o ${name}_tab.csv --remove_all_but filtN_N,filtN_prop_mean,filtNAB_N,filtNAB_prop_mean,TsTv_A,TsTv_B,TsTv_filt,TsTv_filtN,TsTv_filtNAB ;done
#for i in *; do cat $i | sed "s/[0-9]_/0_/g" > $i.mod;done
#for i in *.csv; do mv $i $i.bkp;done
#for i in *.mod; do name=$(echo $i | sed "s/.mod//"); mv $i $name;done

##Output post-processing (problem samples)
#TsTv for problem samples (heter) (outdated)
#for i in ~/dcis/problem_data/run2/DCIS*;do if [[ -d $i ]]; then perl ~/ngcchome/dcis/scripts_vcf/tstv_pretabulate_heter.pl $i/tstv.csv $i/tstv.tomerge; name=$(basename $i);maindir=$(echo $i | sed "s/\/[^\/]*$//"); cat ${maindir}/${name}.csv | awk 'BEGIN{OFS=",";FS=","}{out="";for (i=2;i<=NF;++i){out=out $i OFS};out=substr(out,0,length(out)-1);print out}' > $i/csvtomerge.temp ;perl ~/ngcchome/dcis/scripts_vcf/merge_csv.pl $i/csvtomerge.temp $i/tstv.tomerge $i/${name}_withtstv.temp;rm -f $i/csvtomerge.temp; cat $i/${name}_withtstv.temp | awk "BEGIN{OFS=\",\"}{if (NR==1) {print \"Sample\",\$0} else {print \"$name\",\$0}}" > ${maindir}/${name}_withtstv.csv; rm -f $i/${name}_withtstv.temp ;fi;done

#Heter (all) (outdated)
#sbatch -p private <(echo -e '#!/bin/bash'"\n./HeterAnalyzer_loop.sh /home/dmalload/dcis/problem_data/rerun_samples /home/dmalload/dcis/problem_data/rerun_samples/rerun.list")

#Indexation of bam files
for i in *.bam; do sbatch -p private <(echo -e '#!/bin/bash'"\nmodule load samtools/1.3.1\n samtools index $i");done

#BAM files merging
while read name a b; do sbatch -p private <(echo -e '#!/bin/bash'"\nmodule load bamtools/2.4.0\n bamtools merge -in $a -in $b -out ${name}.bam");done<tomerge.txt

#ychrom
########

#Removes the y chromosome from the human genome
perl -p0e 's/>Y.*?>/>/s' GRCh37-lite.fa > GRCh37-lite_noY.fa

#Fasta index
samtools faidx /home/dmalload/temp_storage/GRCh37-lite_noY.fa

#BWA indexes (needs quite a lot of RAM)
bwa index /home/dmalload/temp_storage/GRCh37-lite.fa
bwa index /home/dmalload/temp_storage/GRCh37-lite_noY.fa

##Checks md5 files in a folder
for i in *.md5; do name=$(echo $i | sed "s/.md5//");echo $(cat $i) " $name" | md5sum -c -;done > out.md5

#Fasta dictionary
module load picard/2.3.0
picard CreateSequenceDictionary R=~/my_storage/GRCh37-lite.fa O=~/my_storage/GRCh37-lite.dict
