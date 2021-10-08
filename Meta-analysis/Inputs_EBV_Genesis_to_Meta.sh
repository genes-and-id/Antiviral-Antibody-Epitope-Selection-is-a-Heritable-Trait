#!/bin/bash
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G
#$ -N Inputs_EBV_Genesis_to_Meta.sh
#$ -o $HOME/logs/EBV_PHENOS/Inputs_EBV_Genesis_to_Meta.out
#$ -e $HOME/logs/EBV_PHENOS/Inputs_EBV_Genesis_to_Meta.err
#$ -t 1-24
#$ -tc 5
echo "**** Job starts ****"
date

## R
module load htslib/1.9
module load R

R CMD BATCH Inputs_EBV_Genesis_to_Meta.R

echo "**** Job ends ****"
date
# job 