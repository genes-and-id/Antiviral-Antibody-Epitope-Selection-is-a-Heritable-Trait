#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G
#$ -N Manual_Manhattan_AFR_Genesis_EBV_Imputed.sh
#$ -o $HOME/logs/2019-05-21_VRC/Manual_Manhattan_AFR_Genesis_EBV_Imputed.out
#$ -e $HOME/logs/2019-05-21_VRC/Manual_Manhattan_AFR_Genesis_EBV_Imputed.err

echo "**** Job starts ****"
date

## R
module load R

R CMD BATCH Manual_Manhattan_AFR_Genesis_EBV_Imputed.R

echo "**** Job ends ****"
date