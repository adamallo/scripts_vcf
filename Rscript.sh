#!/bin/bash

module unload gcc/4.9.2
module load r/3.5.2

Rscript $@
