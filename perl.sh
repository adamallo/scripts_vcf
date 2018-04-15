#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu 7500

###The memory is generating problems

module load perl/5.22.1
perl $@
