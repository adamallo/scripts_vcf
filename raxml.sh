#!/bin/bash

name=$1

#raxmlHPC-AVX2 -f d -m ASC_GTRGAMMA -n ${name}_ml -s $name --asc-corr=lewis -p 2020 -N 20
#raxmlHPC-AVX2 -f d -m ASC_GTRGAMMA -n ${name}_bs -s $name --asc-corr=lewis -p 2020 -b 2020 -N 100
#raxmlHPC-AVX2 -f b -m ASC_GTRGAMMA -n ${name}_bs_tree -z RAxML_bootstrap.${name}_bs -t RAxML_bestTree.${name}_ml
raxmlHPC-AVX2 -f d -m GTRGAMMA -n ${name}_ml -s $name -p 2020 -N 20
raxmlHPC-AVX2 -f d -m GTRGAMMA -n ${name}_bs -s $name -p 2020 -b 2020 -N 100
raxmlHPC-AVX2 -f b -m GTRGAMMA -n ${name}_bs_tree -z RAxML_bootstrap.${name}_bs -t RAxML_bestTree.${name}_ml
