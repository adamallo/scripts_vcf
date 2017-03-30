#!/bin/bash
module load gatk/3.5.0

DATA=/home/dmalload/dcis/problem_data/new_data
GENOME=/home/dmalload/my_storage/GRCh37-lite.fa

cd $DATA

files=("A" "B")
filters=("N" "NAB")
while read ID N A B
do 
    for filter in "${filters[@]}"
    do
        for problem in "${files[@]}"
        do
            for vcffile in $ID/${problem}filt${filter}#*_different.vcf
            do
                for comp in "${files[@]}"
                do
                    if [[ $comp != $problem ]]
                    then
                        echo "Geting the number of alternatives for $comp from $vcffile"
                        vcf2bed --deletions < $vcffile > $ID/${problem}_${filter}_deletions.bed
                        vcf2bed --snvs < $vcffile > $ID/${problem}_${filter}_snvs.bed
                        bedops --everything $ID/${problem}_${filter}_{deletions,snvs}.bed > $ID/${problem}_${filter}.bed
                        if [[ ! -f $ID/${comp}tocompareto${problem}_${filter}.vcf ]]
                        then
                            java -Xms512m -Xmx6G -jar $GATKJAR -T UnifiedGenotyper -R $GENOME -I ${!comp} -o $ID/${comp}tocompareto${problem}_${filter}.vcf --intervals $ID/${problem}_${filter}.bed --output_mode EMIT_ALL_SITES > $ID/${comp}tocompareto${problem}_${filter}.log 2>&1
                        fi
                        cat $ID/${comp}tocompareto${problem}_${filter}.vcf | sed "/^#/d" | perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(",",$F[9]);$reads[1]=="" and $reads[1]=0;print join(",",@F[0,1,3,4],@reads)' > $ID/${problem}_${filter}_counts.csv
                        #save in tabular format
                    fi
                    #else
                       #echo "I cannot compare $problem against $comp"
                    #fi
                done
                
            done
        done
    done
done < listfiles.txt
