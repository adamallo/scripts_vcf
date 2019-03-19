#!/bin/bash

usage="$0 outputname\n"
sufix=".csv"


if [[ $# -ne 1 ]]
then
    echo -e $usage
    exit 1
fi

dir=$(dirname $1)
name=$(basename $1 | sed "s/^\([^.]*\).*$/\1/")
pref=$dir/$name
output="$pref$sufix"
withheader=0

#cat outs
for i in $pref*$sufix
do
    if [[ $withheader -eq 0 ]]
    then
        cat $i > $output
        withheader=1
    else
        tail -n+2 $i >> $output
    fi
    rm -f $i
done

for i in $pref/$name*/
do
    cat $i/vcfdict* | sort | uniq > $i/vcfdict.csv
    cat $i/listdict* | sort | uniq > $i/listdict.csv
    rm -f $i/vcfdict.*.*
    rm -f $i/listdict.*.*
done
