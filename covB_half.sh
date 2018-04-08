#!/bin/bash

usage="$0 vcf outputvcf bam \n Environment variable GATKJAR should point to GATK\'s jar executable"

if [[ $# -ne 3 ]] || [[ ! -f $1 ]] || [[ ! -f $3 ]]
then
    echo -e $usage
    exit 1
fi

vcf=$1
out=$2
bam=$3

module load gatk/3.5.0

GENOME=/home/dmalload/my_storage/GRCh37-lite.fa

dothething ()
{
    bam1=$1
    vcf2=$2
    out=$3

    if [[ ! -f $out ]]
    then
        name_vcf2=$(echo $2 | sed "s/.vcf//g")
        name_out=$(echo $out | sed "s/.tsv//g")
        vcf2bed --deletions < $vcf2 > ${name_vcf2}_deletions.bed
        vcf2bed --snvs < $vcf2 > ${name_vcf2}_snvs.bed
        bedops --everything ${name_vcf2}_{deletions,snvs}.bed | awk 'BEGIN{OFS="\t"}{print($1,$2,$3)}' > ${name_vcf2}.bed
        java -Xms512m -Xmx6G -jar $GATKJAR -T UnifiedGenotyper -R $GENOME -I $bam1 -o "$name_out.vcf" --intervals ${name_vcf2}.bed --output_mode EMIT_ALL_SITES > "$name_out.log" 2>&1
        cat "$name_out.vcf" | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[1]=="" and $reads[1]=0;if($reads[0] eq "./."){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(",",@reads)};print join("\t",@F[0,1,3,4],$readsref,$readsout)' > $out
    else
        echo "The file $out is already present and will be reused"
    fi
}

dothething $bam $vcf $out
