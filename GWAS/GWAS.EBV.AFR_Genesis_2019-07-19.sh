#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -N GWAS.EBV.AFR_Genesis_2019-07-19.sh
#$ -o $HOME/logs/VRC_Entero/GWAS.EBV.AFR_Genesis_2019-07-19.out
#$ -e $HOME/logs/VRC_Entero/GWAS.EBV.AFR_Genesis_2019-07-19.err
#$ -t 1-57
#$ -tc 8
echo "**** Job starts ****"
date

# THIS CODE LUNCH THE R-ARRAY JOB : GWAS.EBV.AFR_Genesis_2019-07-19.R

## R
module load htslib/1.9
module load R

R CMD BATCH GWAS.EBV.AFR_Genesis_2019-07-19.R

echo "**** Job ends ****"
date
# job 