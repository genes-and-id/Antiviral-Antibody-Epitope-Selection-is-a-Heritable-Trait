#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G
#$ -N EBV_Inputs_HLA_mantra.dat.sh
#$ -o $HOME/logs/2019-05-21_VRC/EBV_Inputs_HLA_mantra.dat.out
#$ -e $HOME/logs/2019-05-21_VRC/EBV_Inputs_HLA_mantra.dat.err
#$ -t 1-5
#$ -tc 5
#$ -m e
echo "**** Job starts ****"
date

# THIS CODE LUNCH THE R-ARRAY JOB : EBV_Inputs_HLA_mantra.dat.R

## R
module load R

R CMD BATCH EBV_Inputs_HLA_mantra.dat.R

echo "**** Job ends ****"
date