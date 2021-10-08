#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G
#$ -N EBV_White_Genesis_SuppTable.sh
#$ -o $HOME/logs/2019-05-21_VRC/EBV_EUR_Genesis_SuppTable.out
#$ -e $HOME/logs/2019-05-21_VRC/EBV_EUR_Genesis_SuppTable.err
echo "**** Job starts ****"
date

# THESE GWAS COME FROM THE Deep_Immune_HumanOmni5 PROJECT 
# GENERATED FROM SCRIT: 

## R
module load R

R CMD BATCH EBV_EUR_Genesis_SuppTable.R

echo "**** Job ends ****"
date
# job 