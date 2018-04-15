#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu 7500

###The memory is generating problems

module load parallel/20140822

n_cores=$1

#I have to get the number of cores from the arguments in $@ instead of from SLURM_JOB_CPUS_PER_NODE since if I increase the required memory I get more nodes that the ones I ask for!!!!!

#files=$(ls *.vcf | egrep '^.*filt.*')
files=$(awk 'BEGIN{FS=","}{if($1~/^.*filt.*/){print $2}}' vcfdict.csv)
parallel --delay "0.2" -j $n_cores "echo \"Annotating {1}\"; $SCRIPTSVCF_DIR/annovar.sh {1}" ::: $files
