#!/bin/bash
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G
#$ -N EBV_Meta.Analysis_META_AllGroups.sh
#$ -o $HOME/logs/EBV_Meta.Analysis_META_AllGroups.out
#$ -e $HOME/logs/EBV_Meta.Analysis_META_AllGroups.err
#$ -t 1-24
#$ -tc 6
#$ -m e
echo "**** Job starts ****"
date

cd /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/

## locate files and ID
EURS=`awk '{print $2}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
TWINS=`awk '{print $3}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
AFR=`awk '{print $4}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
LEURs=`awk '{print $5}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
Ltwins=`awk '{print $6}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
Lafr=`awk '{print $7}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
SEURs=`awk '{print $8}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
Stwins=`awk '{print $9}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
Safr=`awk '{print $10}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}"`
ID=$(awk '{print $1}' /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/SampleManifest_METAL_EBV_3groups.txt | awk "NR==${SGE_TASK_ID}")

## Create a temporary directory to do the work 
mkdir -p $TMPDIR/$ID
## cd into the temporary directory 
cd  $TMPDIR/$ID

## RUN META-ANALYSIS USING METHOD 1 'inverse-variance method based on a fixed-effects model'
# SPLIT EACH INPUT FILE INTO CHROMOSOME;
# WE HAVE AUTOSOMAL CHROMOSOMES (1 to 22)
# THEN RUN META ON EACH GROUP OF SPLITTED FILES ;
# COPY THE RESULTS TO THE META-ANALYSIS DIRECTORY 
for i in {1..22};do 
	zcat $EURS | awk -v chr=$i '{if (NR==1 || $1==chr) print}' | gzip -c > $ID.chr$i.EURs.txt.gz;
	zcat $TWINS | awk -v chr=$i '{if (NR==1 || $1==chr) print}' | gzip -c > $ID.chr$i.twins.txt.gz;
	zcat $AFR | awk -v chr=$i '{if (NR==1 || $1==chr) print}' | gzip -c > $ID.chr$i.afr.txt.gz; 
	~/toolsCV/META/meta --method 1 --threshold 0.3 --lambda $LEURs $Ltwins $Lafr --cohort $ID.chr$i.EURs.txt.gz $ID.chr$i.twins.txt.gz $ID.chr$i.afr.txt.gz --output $ID.chr$i.txt;
	cp $ID.chr$i.txt /dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/META_ANALYSIS/ ;
	done 

## CLEAN 
rm -rf $TMPDIR/$ID



echo "**** Job ends ****"
date