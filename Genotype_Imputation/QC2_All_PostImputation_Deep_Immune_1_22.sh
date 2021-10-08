#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=30G
#$ -N QC2_All_PostImputation_Deep_Immune_1_22.sh
#$ -o $HOME/logs/VRC_Entero/QC2_All_PostImputation_Deep_Immune_1_22.out
#$ -e $HOME/logs/VRC_Entero/QC2_All_PostImputation_Deep_Immune_1_22.err
echo "**** Job starts ****"
date

module load vcftools/0.1.14
module load htslib/1.9
module load R

## RUN 
R CMD BATCH QC2_All_PostImputation_Deep_Immune_1_22.R

# done
echo "**** Job ends ****"
date