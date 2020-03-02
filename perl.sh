#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH --mem-per-cpu 7500

###The memory is generating problems
module load perl/5.26.0
module load htslib/1.2.1
eval "$(perl -I/home/afortun2/perl5/lib/perl5 -Mlocal::lib=/home/afortun2/perl5)" #My own perl library
#perl 5.26.0 loaded perl libraries added to @INC
perl $@
