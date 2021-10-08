#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G
#$ -N EBV_dmatcal_MANTRA_HLA.sh
#$ -o $HOME/logs/EBV_MANTRA_HLA.out
#$ -e $HOME/logs/EBV_MANTRA_HLA.err
#$ -m e
echo "**** Job starts ****"
date

cd ~/scratch

## Create a temporary directory to do the work 
mkdir -p TMP
## cd into the temporary directory 
cd  TMP

## COPY FILES 
cp /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/mantra.in ./
cp /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/33024_HLA_mantra.dat ./

## Change name for this step 
mv 33024_HLA_mantra.dat mantra.dat 

## run dmatcal 
~/toolsCV/MANTRA_software/dmatcal
## copy file
cp dmat.out /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/
### CLEAN 
rm -rf ~/scratch/TMP

## DONE 
echo "**** Job ends ****"
date