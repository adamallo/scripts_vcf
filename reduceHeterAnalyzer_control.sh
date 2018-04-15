#!/bin/bash

usage="$0 outputname\n"
sufix=".csv"


if [[ $# -ne 1 ]]
then
    echo -e $usage
    exit 1
fi

pref=$(echo $1 | sed "s/^\([^.]*\).*$/\1/")
output="$pref$sufix"
withheader=0

#cat outs
for i in $pref*
do
    if [[ $withheader -eq 0 ]]
    then
        cat $i > $output
        withheader=1
    else
        tail -n+2 $i >> $output
    fi
done

cd $pref

cat vcfdict* | sort | uniq > vcfdict.csv
cat listdict* | sort | uniq > listdict.csv
