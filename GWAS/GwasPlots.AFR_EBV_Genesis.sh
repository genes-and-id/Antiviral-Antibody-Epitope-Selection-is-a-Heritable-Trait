#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G
#$ -N GwasPlots.AFR_EBV_Genesis.sh
#$ -o $HOME/logs/2019-05-21_VRC/GwasPlots.AFR_EBV_Genesis.out
#$ -e $HOME/logs/2019-05-21_VRC/GwasPlots.AFR_EBV_Genesis.err
#$ -t 1-27
#$ -tc 8
echo "**** Job starts ****"
date

# THIS CODE LUNCH THE R-ARRAY JOB : GwasPlots.AFR_EBV_Genesis.R
# TO GENERATE MANHATTAN PLOTS AND Q-Q PLOTS FROM THE MULTIPLE GWAS
# THESE GWAS COME FROM THE Deep_Immune_HumanOmni5 PROJECT 

## Generate ID
ls -d /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis/*.assoc > ~/scratch/EBV.AFR.list
ID=`awk '{print $1}' ~/scratch/EBV.AFR.list | awk "NR==${SGE_TASK_ID}"`

## R
module load R

R CMD BATCH GwasPlots.AFR_EBV_Genesis.R

echo "**** Job ends ****"
date