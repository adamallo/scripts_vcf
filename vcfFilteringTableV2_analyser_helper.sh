#!/bin/bash
#
cd $PBS_O_WORKDIR
module load perl/5.22.1

perl $@ --n_cores $PBS_NUM_PPN
