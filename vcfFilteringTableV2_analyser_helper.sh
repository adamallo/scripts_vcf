#!/bin/bash
#
#SBATCH -p private
#SBATCH -N 1
#SBATCH -n 1

##module load perl/5.22.1 Current module is broken in the cluster
module load parallel/20140822

perl $@ --n_cores $SLURM_JOB_CPUS_PER_NODE

ordir=${@[10]}

if [[ -x $ordir/annovar.sh ]]
then

	for i in *filterN*.vcf
	do
		sem -j+0 echo "Annotating $i" ";" $ordir/annovar.sh $i
	done
	sem --wait
else
	echo "Error finding the executable annovar.sh"
	exit
fi
