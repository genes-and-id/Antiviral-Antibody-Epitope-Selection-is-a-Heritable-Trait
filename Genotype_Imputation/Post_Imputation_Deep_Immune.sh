#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=50G
#$ -N Post_Imputation_Deep_Immune.sh
#$ -o $HOME/logs/Post_Imputation_Deep_Immune.out
#$ -e $HOME/logs/Post_Imputation_Deep_Immune.err
echo "**** Job starts ****"
date


module load bcftools/1.2
module load vcftools/0.1.14
module load htslib/1.9 


## FILES
NFILE=Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.vcf.gz
TXT=/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.txt
VARS=/users/cvalenci/toolsCV/dbSNP_reference/All_20170710.GRCh37p13.vcf.gz

cd /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/OUTPUT_MIS/VCF

## make sinlge VCF 
bcftools concat -f /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/List_vcf.txt | bgzip -c > temp.vcf.gz

## index 
tabix -p vcf temp.vcf.gz

## Annotate variant ID
bcftools annotate -a $VARS -c ID --output $NFILE --output-type z temp.vcf.gz

# Index final VCF 
tabix -p vcf $NFILE 

## Get numbers  
zcat $NFILE | grep -v '^#' | awk '{OFS="\t"; split($8,x,';'); print $1,$2,$3,$7,x[1],x[2],x[3],x[4]}' > $TXT

## clean
rm /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/OUTPUT_MIS/VCF/temp.vcf.gz*

echo "**** Job ends ****"
date