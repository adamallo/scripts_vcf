#!/bin/bash
#module load bcftools/1.4.0

#echo "sample;#AN;#BN;#comun;sim" > results.noscript.csv
#for i in *
#do 
    i=$1
    if [[ -d $i ]]
    then 
        
        #Filtering A and b considering only PASS and allelebias

        java -Xmx2G -jar $SNPSIFT_DIR/SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS') | (FILTER = 'alleleBias')" $i/A.vcf > $i/A_filt.vcf
        java -Xmx2G -jar $SNPSIFT_DIR/SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS') | (FILTER = 'alleleBias')" $i/B.vcf > $i/B_filt.vcf

        #Intersection of A and B
        bgzip $i/A_filt.vcf; bgzip $i/B_filt.vcf;tabix $i/A_filt.vcf.gz;tabix $i/B_filt.vcf.gz;bcftools isec $i/A_filt.vcf.gz $i/B_filt.vcf.gz -p $i/bcfcomp

        #Intersection of N and CAB (all variants, CAB is only PASS and allelebias, we want to be lenient with N)
        bgzip $i/N.vcf; bgzip $i/bcfcomp/0002.vcf; tabix $i/N.vcf.gz; tabix $i/bcfcomp/0002.vcf.gz;bcftools isec $i/bcfcomp/0002.vcf.gz $i/N.vcf.gz -p $i/bcfcomp/CAP
        
        #Intersection of N and PA
        bgzip $i/bcfcomp/0000.vcf; tabix $i/bcfcomp/0000.vcf.gz;bcftools isec $i/bcfcomp/0000.vcf.gz $i/N.vcf.gz -p $i/bcfcomp/PA
        
        #Intersection of N and PB
        bgzip $i/bcfcomp/0001.vcf; tabix $i/bcfcomp/0001.vcf.gz;bcftools isec $i/bcfcomp/0001.vcf.gz $i/N.vcf.gz -p $i/bcfcomp/PB

        #Calculation of %prop and Nvariants
        na=$(cat $i/bcfcomp/PA/0000.vcf | grep -v "#" | wc -l)
        nb=$(cat $i/bcfcomp/PB/0000.vcf | grep -v "#" | wc -l)
        ncom=$(cat $i/bcfcomp/CAP/0000.vcf | grep -v "#" | wc -l)
        tot=$(($na+$nb+$ncom))
        prop=$(echo "scale=3;$ncom/$tot" | bc)

        echo "$i;$na;$nb;$ncom;$prop"   
    fi
#done >> results.noscript.csv
#done

