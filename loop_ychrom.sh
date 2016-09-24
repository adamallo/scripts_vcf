#!/bin/bash
#SBATCH -t 1-00:00

refgenome="/home/dmalload/temp_storage/GRCh37-lite.fa"
refgenomenoy="/home/dmalload/temp_storage/GRCh37-lite_noY.fa"
bedfile="/home/dmalload/ngcchome/xgen-exome-research-panel-targets.bed"
usage="$0 directory torun_file min_mapq max_mapq_reject"

if [[ $# -ne 4 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]
then
    echo $usage
    exit
else
    dir=$(readlink -f $1)
    torun=$(readlink -f $2)
    min_mapq=$3
    max_mapq=$4
fi

while read -r output normal a b
do
sbatch -p private --mem=20000 <(echo -e '#!/bin/bash'"\n perl $SCRIPTSVCF_DIR/ychrom.pl $normal $refgenome $refgenomenoy $bedfile $min_mapq $max_mapq $dir ${output}_N > $dir/${output}_N.out")
sbatch -p private --mem=20000 <(echo -e '#!/bin/bash'"\n perl $SCRIPTSVCF_DIR/ychrom.pl $a $refgenome $refgenomenoy $bedfile $min_mapq $max_mapq $dir ${output}_A > $dir/${output}_A.out")
sbatch -p private --mem=20000 <(echo -e '#!/bin/bash'"\n perl $SCRIPTSVCF_DIR/ychrom.pl $b $refgenome $refgenomenoy $bedfile $min_mapq $max_mapq $dir ${output}_B > $dir/${output}_B.out")
sleep 30
done < $torun
