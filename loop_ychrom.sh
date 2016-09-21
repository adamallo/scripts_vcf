#!/bin/bash
#SBATCH -t 4-00:00

refgenome="/home/dmalload/temp_storage/GRCh37-lite.fa"
bedfile="/home/dmalload/ngcchome/xgen-exome-research-panel-targets.bed"
usage="$0 directory torun_file min_mapq"

if [[ $# -ne 3 ]] || [[ ! -d $1 ]] || [[ ! -f $2 ]]
then
    echo $usage
    exit
else
    dir=$(readlink -f $1)
    torun=$(readlink -f $2)
    min_mapq=$3
fi

while read -r output normal a b
do
sbatch -p private --mem=20000 <(echo -e '#!/bin/bash'"\n perl $SCRIPTSVCF_DIR/ychrom.pl $normal $refgenome $bedfile $min_mapq > $dir/${output}_N.out")
sbatch -p private --mem=20000 <(echo -e '#!/bin/bash'"\n perl $SCRIPTSVCF_DIR/ychrom.pl $a $refgenome $bedfile $min_mapq > $dir/${output}_A.out")
sbatch -p private --mem=20000 <(echo -e '#!/bin/bash'"\n perl $SCRIPTSVCF_DIR/ychrom.pl $b $refgenome $bedfile $min_mapq > $dir/${output}_B.out")
done < $torun
