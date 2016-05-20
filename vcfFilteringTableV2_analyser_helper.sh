#!/bin/bash
#
cd $PBS_O_WORKDIR
module load perl/5.22.1

perl $@ --n_cores $PBS_NUM_PPN

ordir=${@[10]}

if [[ -x $ordir/annovar.sh ]]
then

	for i in *filterN*.vcf
	do
		sem -j+0 $ordir/annovar.sh $i
	done
	sem --wait
else
	echo "Error finding the executable annovar.sh"
	exit
fi
