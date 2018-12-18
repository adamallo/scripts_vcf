#!/bin/bash

usage="$0 outputname\n"
sufix=".csv"


if [[ $# -ne 1 ]]
then
    echo -e $usage
    exit 1
fi

dir=$(dirname $1)
pref=$(basename $1 | sed "s/^\([^.]*\).*$/\1/")
pref=$dir/$pref
output="$pref$sufix"
withheader=0

#cat outs
for i in $pref.*$sufix
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

cd $pref

cat vcfdict* | sort | uniq > vcfdict.csv
cat listdict* | sort | uniq > listdict.csv
rm -f vcfdict.*.*
rm -f listdict.*.*
