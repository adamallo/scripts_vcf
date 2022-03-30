#!/bin/bash
module load bedtools2/2.24.0

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


# FOREACH FOLDER (PATIENT)
for patientLog in $mainDir/*.out
do
    patient=$(basename $patientLog | sed "s/.out//")
    workDir="$mainDir/$patient"
    dataFileList="$workDir/vcfdict.csv"

    sortedListOfGenes=$(echo $(perl -lane '{print $F[3]}' $bedFile | sort | uniq))
    nGenes=( $sortedListOfGenes )
    nGenes=${#nGenes[@]}


    #Initial
        #This will require calculating the union of two vcf files
    
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

    printFun "$sortedListOfGenes" "$resultF" $patient "F"
    printFun "$sortedListOfGenes" "$resultFS" $patient "FS"
    printFun "$sortedListOfGenes" "$resultFSB" $patient "FSB"
    printFun "$sortedListOfGenes" "$resultFSBN" $patient "FSBN"
    printFun "$sortedListOfGenes" "$resultFSBNP" $patient "FSBNP"
     
done

