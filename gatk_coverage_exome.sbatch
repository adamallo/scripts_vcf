#!/bin/bash
#SBATCH -n1
#SBATCH -c1
#SBATCH --mem-per-cpu 10000

module load gatk/3.5.0
name=$1
java -Xms512m -Xmx6G -jar /packages/6x/gatk/3.5.0/GenomeAnalysisTK.jar -T DepthOfCoverage -R $HUMAN_GENOME -L /home/dmalload/dcis/beds/xgen-exome-research-panel-targets.bed -I $2 -o ${name}_coverage --summaryCoverageThreshold 100 --summaryCoverageThreshold 80 --summaryCoverageThreshold 60 --summaryCoverageThreshold 50 --summaryCoverageThreshold 40 --summaryCoverageThreshold 30 --summaryCoverageThreshold 20 --summaryCoverageThreshold 15 --summaryCoverageThreshold 10
