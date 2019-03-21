#!/bin/bash

usage="$0 bed outputvcf bam \n Environment variable GATKJAR should point to GATK\'s jar executable"

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ ! -f $3 ]]
then
    echo -e $usage
    exit 1
fi

bed=$1
out=$2
bam=$3

module load gatk/3.5.0

GENOME=/home/dmalload/my_storage/GRCh37-lite.fa

dothething ()
{
    bam1=$1
    bed=$2
    out=$3

    if [[ ! -f $out ]]
    then
        name_out=$(echo $out | sed "s/.tsv//g")
        java -Xms512m -Xmx6G -jar $GATKJAR -T UnifiedGenotyper -R $GENOME -I $bam1 -o "$name_out.vcf" --intervals $bed --output_mode EMIT_ALL_SITES -glm BOTH > "$name_out.log" 2>&1
        cat "$name_out.vcf" | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[1]=="" and $reads[1]=0;if($reads[0] eq "./."){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(",",@reads)};print join("\t",@F[0,1,3,4],$readsref,$readsout)' > $out
    else
        echo "The file $out is already present and will be reused"
    fi
}

dothething $bam $bed $out
