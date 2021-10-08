#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G
#$ -N EBV_MANTRA_HLA.sh
#$ -o $HOME/logs/EBV_MANTRA_HLA.out
#$ -e $HOME/logs/EBV_MANTRA_HLA.err
#$ -t 1-5
#$ -tc 5
#$ -m e
echo "**** Job starts ****"
date

## MANTRA ANALYSIS:
## 3. RUN MANTRA ANALYSIS IN HLA
## HOW TO RUN: qsub EBV_Inputs_HLA_mantra.dat.sh 

## Locate files and IDs
FILE=`awk '{print $1}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/SampleManifest_EBV_MANTRA_HLA.txt | awk "NR==${SGE_TASK_ID}"`
NAME=`awk '{print $3}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/SampleManifest_EBV_MANTRA_HLA.txt | awk "NR==${SGE_TASK_ID}"`
ID=$(awk '{print $2}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/SampleManifest_EBV_MANTRA_HLA.txt | awk "NR==${SGE_TASK_ID}")

## define directories 
mkdir -p /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/HLA_Results/
#
cd /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/HLA_Results/

## Create a temporary directory to do the work 
mkdir -p $ID
pwd
## cd into the temporary directory 
cd $ID

## copy files: Input.dat, mantra.in, and mantra.out (1 for all as Fst should be the same)
cp $FILE ./
cp /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/mantra.in ./
cp /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/dmat.out ./
cp /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/fname.in ./

## CHANGE NAME INPUT DAT (SERIOUSLY??)
mv $NAME mantra.dat 

## RUN MANTRA 
~/toolsCV/MANTRA_software/mantra.v1 <fname.in 


## DONE 
echo "**** Job ends ****"
date