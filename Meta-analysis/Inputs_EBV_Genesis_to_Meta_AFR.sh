#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -N Inputs_EBV_Genesis_to_Meta_AFR.sh
#$ -o $HOME/logs/EBV_PHENOS/Inputs_EBV_Genesis_to_Meta_AFR.out
#$ -e $HOME/logs/EBV_PHENOS/Inputs_EBV_Genesis_to_Meta_AFR.err
#$ -t 1-41
#$ -tc 3
echo "**** Job starts ****"
date


## R
module load htslib/1.9
module load R

R CMD BATCH Inputs_EBV_Genesis_to_Meta_AFR.R

echo "**** Job ends ****"
date
# job 