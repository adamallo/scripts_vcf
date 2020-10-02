#!/bin/bash
###-t 4-00:00
###--mem-per-cpu 7500

###The memory is generating problems
module load perl/5.26.0
eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)" #My own perl library
#perl 5.26.0 loaded perl libraries added to @INC
perl $@
