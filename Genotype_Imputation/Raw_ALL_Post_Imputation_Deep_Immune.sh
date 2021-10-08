#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=50G
#$ -N Raw_ALL_Post_Imputation_Deep_Immune.sh
#$ -o $HOME/logs/Raw_All_Post_Imputation_Deep_Immune.out
#$ -e $HOME/logs/Raw_All_Post_Imputation_Deep_Immune.err
echo "**** Job starts ****"
date

## FIRST STEP IN QC AFTER IMPUTATION FILE THAT INCLUDES ALL VRC PARTICIPANTS
## THEN CALCULATE ALLELE FREQ
## AND EXTRACT R-SQUARE 

module load vcftools/0.1.14
module load htslib/1.9 

## FILES 
VCF=/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/OUTPUT_MIS/VCF/Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.vcf.gz
FQ=/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/rawALL_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe
TXT=/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/rawALL_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.txt

## get allele freq 
vcftools --gzvcf $VCF --freq --out /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/rawALL_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe 

## get R2 , Note split always with double quote " "
zcat $VCF | grep -v '^#' | awk '{OFS="\t"; split($8,x,";"); print $1,$2,$3,$7,x[1],x[2],x[3],x[4]}' > $TXT

echo "**** Job ends ****"
date