#!/bin/bash
module load bcftools/1.12.0

mainDir="/scratch/dmalload/dcis/aim4Genes/data"
bedFile="/scratch/dmalload/dcis/aim4Genes/finalGenes.bed"

##FUNCTIONS

# Uses a pattern and a vcfdict file and returns a space-separated string of 2 files or -1 on error
function getVCFs {
    pattern=$1
    dataFileList=$2
    workDir=$(dirname $2)
    
    if [[ ! -f $dataFileList ]]
    then
        echo -1
        break
    fi

    vcfFiles=( $(grep -P -e $pattern $dataFileList | perl -F, -slane '{print "$workDir"."/".$F[1]}' -- -workDir="$workDir") )

    if [[ ${#vcfFiles[@]} -ne 2 ]]
    then
        echo -1
        break
    fi

    echo "${vcfFiles[@]}"
}

# Uses a string of 2 space-separated vcf files and a bedfile, and returns a string of space-separated values with the number of vcf elements in each of the interval (genes) in the bed file
function getCounts {
    bedFile=$1
    vcfFiles=( $2 )

    #echo -e "$(bedtools intersect -a $bedFile -b ${vcfFiles[@]} -wb | perl -lane '{print $F[3]}')""\n""$(perl -lane '{print $F[3]}' $bedFile)" | sort | uniq -c | perl -lane 'BEGIN{$,="\t"}{print($F[1],$F[0]-1)}' 1>&2 #Old version with two columns
    
    #To count the variants for each gene keeping all of them, I am concatenating the output of the bedtools intersect with the list of genes, sorting and counting the number of each of them, and then substracting one
    echo $(echo -e "$(bedtools intersect -a $bedFile -b ${vcfFiles[@]} -wb | perl -lane '{print $F[3]}')""\n""$(perl -lane '{print $F[3]}' $bedFile)" | sort | uniq -c | perl -lane '$F[1] =~ /^\s*$/ or {print($F[0]-1)}')
}

# Uses a string of 2 space-separated vcf files and a bedfile, and returns a string of space-separated values with the number of vcf elements in each of the interval (genes) in the bed file
# In this case, the VCF files need to be filtered for PASS || alleleBias and may contain overlapping entries
function getCountsUnprocessed {
    bedFile=$1
    vcfFiles=( $2 )

    #To count the variants for each gene keeping all of them, I am concatenating the output of the bedtools intersect with the list of genes, sorting and counting the number of each of them, and then substracting one
    echo $(echo -e "$(echo "$(bcftools view -i 'FILTER="PASS" || FILTER="alleleBias"' ${vcfFiles[0]} 2>/dev/null)""\n""$(bcftools view -i 'FILTER="PASS" || FILTER="alleleBias"' ${vcfFiles[1]} 2>/dev/null)" | perl -lane 'BEGIN{$,="\t"}{/^#/ or print($F[0],$F[1]-1,$F[1])}' | sort -k1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -a $bedFile -b stdin -wb | perl -lane '{print $F[3]}')""\n""$(perl -lane '{print $F[3]}' $bedFile)" | sort | uniq -c | perl -lane '$F[1] =~ /^\s*$/ or {print($F[0]-1)}')
}

# Uses a vcf file and a bedfile, and returns a string of space-separated values with the number of vcf elements in each of the interval (genes) in the bed file
function getCountsUnprocessedSingleVCF {
    bedFile=$1
    vcfFile=$2

    #To count the variants for each gene keeping all of them, I am concatenating the output of the bedtools intersect with the list of genes, sorting and counting the number of each of them, and then substracting one
    echo $(echo -e "$(perl -lane 'BEGIN{$,="\t"}{/^#/ or print($F[0],$F[1]-1,$F[1])}' $vcfFile | sort -k1,1 -k2,2n | bedtools intersect -a $bedFile -b stdin -wb | perl -lane '{print $F[3]}')""\n""$(perl -lane '{print $F[3]}' $bedFile)" | sort | uniq -c | perl -lane '$F[1] =~ /^\s*$/ or {print($F[0]-1)}')
}

# Uses a vcf file and a bedfile, and returns a string of space-separated values with the number of vcf elements in each of the interval (genes) in the bed file
function getCountsFilteredSingleVCF {
    bedFile=$1
    vcfFile=$2

    #To count the variants for each gene keeping all of them, I am concatenating the output of the bedtools intersect with the list of genes, sorting and counting the number of each of them, and then substracting one
    echo $(echo -e "$(bcftools view -i 'FILTER="PASS" || FILTER="alleleBias"' $vcfFile 2>/dev/null | perl -lane 'BEGIN{$,="\t"}{/^#/ or print($F[0],$F[1]-1,$F[1])}' | sort -k1,1 -k2,2n | bedtools intersect -a $bedFile -b stdin -wb | perl -lane '{print $F[3]}')""\n""$(perl -lane '{print $F[3]}' $bedFile)" | sort | uniq -c | perl -lane '$F[1] =~ /^\s*$/ or {print($F[0]-1)}')
}

# Uses a string of 2 space-separated vcf files a bedfile, and a VCF file with variants to discard, and returns a string of space-separated values with the number of vcf elements in each of the interval (genes) in the bed file
# In this case, the VCF files (not the ones to discard) need to be filtered for PASS || alleleBias and may contain overlapping entries
function getCountsUnprocessedDiscardVars {
    bedFile=$1
    vcfFiles=( $2 )
    vcfN=$3

    #To count the variants for each gene keeping all of them, I am concatenating the output of the bedtools intersect with the list of genes, sorting and counting the number of each of them, and then substracting one
    echo $(echo -e "$(echo "$(bcftools view -i 'FILTER="PASS" || FILTER="alleleBias"' ${vcfFiles[0]} 2>/dev/null)""\n""$(bcftools view -i 'FILTER="PASS" || FILTER="alleleBias"' ${vcfFiles[1]} 2>/dev/null)" | perl -lane 'BEGIN{$,="\t"}{/^#/ or print($F[0],$F[1]-1,$F[1])}' | sort -k1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -v -wa -a stdin -b $vcfN | bedtools intersect -a $bedFile -b stdin -wb | perl -lane '{print $F[3]}')""\n""$(perl -lane '{print $F[3]}' $bedFile)" | sort | uniq -c | perl -lane '$F[1] =~ /^\s*$/ or {print($F[0]-1)}')
}
    
# Gets space-separated strings with gene names,and result, and names of the patient and filtering-stage and prints the data in long format
function printFun {
    sortedListOfGenes=$1
    resultString=$2
    patient=$3
    filter=$4
    temp=( $sortedListOfGenes )
    n=${#temp[@]}
    sep=" "

    filterString=$(perl -se 'print(join($sep,($text) x $times))' -- -sep=" " -text=$filter -times=$n)
    patientString=$(perl -se 'print(join($sep,($text) x $times))' -- -sep=" " -text=$patient -times=$n)
    echo -e "$patientString\n$filterString\n$sortedListOfGenes\n$resultString" | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{ii=1;for(i in a){jj=1;for(j in a[ii]){printf"%s"(jj==NR?"\n":FS),a[ii][jj];jj++};ii++}}' FS=" "
}

# MAIN

usageMessage="\nUsage: $0 [-h] [-b bedfile] [-d dir].\n\nDefaults:\n\tbedfile: $bedFile\n\tdir: $mainDir\n"
function usage {
    if [[ $# -eq 1 ]]
    then
        echo -e "\nERROR: $1\n""$usageMessage" 1>&2
    else
        echo -e "$usageMessage" 1>&2
    fi
    
    exit 1
}

#The first : indicates silent checking
#The : after each option indicates that they require an argument
while getopts ":b:d:h" options
do
    case "${options}" in 
        b)
            bedFile=${OPTARG}
            [ -f "$bedFile" ] || usage "bed file not found"
            ;;
        d)
            mainDir=${OPTARG}
            [ -d "$mainDir" ] || usage "input directory not found"
            ;;
        h)
            usage
            ;;
        *)
            usage "wrong input option"
            ;;
    esac
done

#After this, there would be positional arguments available
shift "$((OPTIND-1))"
[ $# -eq 0 ] || usage "this script does not allow positional arguments"

#Let's top playing with getopts and now this is the real program

# FOREACH FOLDER (PATIENT)
for patientLog in $mainDir/*.out
do
    patient=$(basename $patientLog | sed "s/.out//")
    workDir="$mainDir/$patient"
    dataFileList="$workDir/vcfdict.csv"

    sortedListOfGenes=$(echo $(perl -lane '{print $F[3]}' $bedFile | sort | uniq))
    nGenes=( $sortedListOfGenes )
    nGenes=${#nGenes[@]}


    ##Initial ones use unprocessed VCF files and thus are different than the rest

    #Initial
    ########
    pattern="^[A-B]#[^#]*#[^#]*.vcf,"
    vcfFiles=$(getVCFs $pattern $dataFileList $workDir)

    if [[ "$vcfFiles" == -1 ]]
    then
        echo "Skipping the directory/patient ${workDir}. VCF files not found"
        continue
    fi

    resultI=$(getCountsUnprocessed $bedFile "$vcfFiles")
    ########
    
    #Germline
    #########
    vcfN="$workDir/N.vcf"

    if [[ ! -f $vcfN ]]
    then
        echo "Skipping the directory/patient $workDir, N.vcf file not found"
        continue
    fi

    resultG=$(getCountsUnprocessedSingleVCF $bedFile $vcfN)
    #########
    
    #Filtered germline
    ##################
    vcfN="$workDir/N.vcf"
    resultFG=$(getCountsFilteredSingleVCF $bedFile $vcfN)
    ##################
    
    #Initial without somatic
    ########################
    vcfN="$workDir/N.vcf"
    resultIS=$(getCountsUnprocessedDiscardVars $bedFile "$vcfFiles" "$vcfN")
    ########################
    
    #Filtered Variants
    ##################
    pattern="^filtU#.*.vcf"
    vcfFiles=$(getVCFs $pattern $dataFileList $workDir)

    if [[ "$vcfFiles" == -1 ]]
    then
        echo "Skipping the directory/patient $workDir"
        continue
    fi

    resultF=$(getCounts $bedFile "$vcfFiles")
    ##################
    
    #Filtered-Somatic variants
    ##########################
    pattern="^filtNU#.*.vcf"
    vcfFiles=$(getVCFs $pattern $dataFileList $workDir)

    if [[ "$vcfFiles" == -1 ]]
    then
        echo "Skipping the directory/patient $workDir"
        continue
    fi

    resultFS=$(getCounts $bedFile "$vcfFiles")
    ##########################

    #Filtered-Somatic-CovB variants
    ###############################
    pattern="^filtNcovBU#.*.vcf"
    vcfFiles=$(getVCFs $pattern $dataFileList $workDir)

    if [[ "$vcfFiles" == -1 ]]
    then
        echo "Skipping the directory/patient $workDir"
        continue
    fi

    resultFSB=$(getCounts $bedFile "$vcfFiles")
    ###############################
    

    #Filtered-Somatic-CovB-CovN variants
    ####################################
    pattern="^filtcovBNABU#.*.vcf"
    vcfFiles=$(getVCFs $pattern $dataFileList $workDir)

    if [[ "$vcfFiles" == -1 ]]
    then
        echo "Skipping the directory/patient $workDir"
        continue
    fi

    resultFSBN=$(getCounts $bedFile "$vcfFiles")
    ####################################
    
    #Filtered-Somatic-CovB-CovN-PAF variants
    ########################################
    pattern="^filtcovBNABUPAF#.*.vcf"
    vcfFiles=$(getVCFs $pattern $dataFileList $workDir)

    if [[ "$vcfFiles" == -1 ]]
    then
        echo "Skipping the directory/patient $workDir"
        continue
    fi

    resultFSBNP=$(getCounts $bedFile "$vcfFiles")
    ########################################

    #AWK TRANSPOSE ONELINER FROM https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash. For some reason it did not work in my awk and I needed to change some things (using increasing counters because the indexes were not being returned sorted)
    #original
    #echo -e "$sortedListOfGenes\n$resultF\n$resultFS\n$resultFSB\n$resultFSBN\n$resultFSBNP" | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{for(i in a)for(j in a[i])printf"%s"(j==NR?"\n":FS),a[i][j]}' FS=" "
    #Fixed
    #echo -e "$sortedListOfGenes\n$resultF\n$resultFS\n$resultFSB\n$resultFSBN\n$resultFSBNP" | awk '{for(i=1;i<=NF;i++)a[i][NR]=$i}END{ii=1;for(i in a){jj=1;for(j in a[ii]){printf"%s"(jj==NR?"\n":FS),a[ii][jj];jj++};ii++}}' FS=" "
    #Long version

    printFun "$sortedListOfGenes" "$resultI" $patient "I"
    printFun "$sortedListOfGenes" "$resultIS" $patient "IS"
    printFun "$sortedListOfGenes" "$resultG" $patient "G"
    printFun "$sortedListOfGenes" "$resultFG" $patient "FG"
    printFun "$sortedListOfGenes" "$resultF" $patient "F"
    printFun "$sortedListOfGenes" "$resultFS" $patient "FS"
    printFun "$sortedListOfGenes" "$resultFSB" $patient "FSB"
    printFun "$sortedListOfGenes" "$resultFSBN" $patient "FSBN"
    printFun "$sortedListOfGenes" "$resultFSBNP" $patient "FSBNP"
     
done

