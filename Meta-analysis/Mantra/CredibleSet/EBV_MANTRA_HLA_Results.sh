#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G
#$ -N EBV_MANTRA_HLA_Results.sh
#$ -o $HOME/logs/2019-05-21_VRC/EBV_MANTRA_HLA_Results.out
#$ -e $HOME/logs/2019-05-21_VRC/EBV_MANTRA_HLA_Results.err
#$ -m e
echo "**** Job starts ****"
date

## R
module load R

R CMD BATCH EBV_MANTRA_HLA_Results.R

echo "**** Job ends ****"
date