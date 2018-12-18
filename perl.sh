#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu 7500

###The memory is generating problems
#module load perl/5.26.0
#perl 5.26.0 is being loaded in my .basrc and my perl libraries added to @INC
perl $@
