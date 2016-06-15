#!/bin/bash
#
#SBATCH -p private
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 4-00:00

module load perl/5.22.1
module load parallel/20140822

args=($@)
last_arg=$(( $#-1 ))
ordir=${args[$last_arg]}

if [[ $# -ne 1 ]]
then
    echo "Variant analysis"
    perl $@ --n_cores $SLURM_JOB_CPUS_PER_NODE
else
    echo "Skipping the analysis of variants, directly executing annovar"
fi

if [[ -x $ordir/annovar.sh ]]
then

    files=$(ls filt*.vcf)
    parallel --delay "0.2" -j $SLURM_JOB_CPUS_PER_NODE --joblog annovar.log --resume "echo \"Annotating {1}\"; $ordir/annovar.sh {1}" ::: $files

#	for i in filt*.vcf
#	do
#		#sem -j$SLURM_JOB_CPUS_PER_NODE echo "Annotating $i" ";" $ordir/annovar.sh $i
#        echo "Anotating $i\n"
#        $ordir/annovar.sh $i
#	done
#	#sem --wait
else
	echo "Error finding the executable annovar.sh"
	exit
fi
