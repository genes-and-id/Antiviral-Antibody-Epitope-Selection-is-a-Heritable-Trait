# QC2_All_PostImputation_Deep_Immune_1_22.R

# module load vcftools/0.1.14
# module load htslib/1.9
# ml R

# TO RUN: qsub QC2_All_PostImputation_Deep_Immune_1_22.sh

## SECOND STEP IN QC AFTER IMPUTATION FILE THAT INCLUDES ALL VRC PARTICIPANTS
## USES FREQ AND Rsq TO FILTER VARIANTS WITH MAF >= 1% and RSQ >=0.3

## INPUT VCF WAS: Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.vcf.gz
## OUTPUT VCF: All_QC.Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.vcf.gz
## OUTPUT GDS: All_QC.Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.gds
## IS ALSO GENERATES THE PHENOTYPE FILE FOR GENESIS (GWAS)


# TO RUN : qsub PCair_All_FiltDeep_Immune_1_22_ibdmafhwe.sh 

# Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(SeqArray)
library(data.table)
data.table::setDTthreads(7) # From Leo:to avoid using too many cores
options(stringsAsFactors=FALSE)


# Create directory for ANALYSIS
main_dir <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/'
sub_dir <- 'All_Participants'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/All_Participants/'

# LOAD INFO (PCA)
PCs <- read.table('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/ph53214_pca_todos.txt',header=T)
PCs <- PCs[,c(2:12)] # extract top 10 PCs
names(PCs)[1] <- 'scanID'

# Pheno for GWAS
pheno <- fread('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_pheno2019-07-19_HumanOmmni.txt',header=T)
df <- as.data.frame(pheno[,-1])
# Keep only PCs for those in pheno
df <- merge(PCs,df,by.x='scanID',by.y='IID') #
df <- df[complete.cases(df),] # 

## SAVE FILE
save(df,file=paste0(file.path,'EBV_pheno2019-07-19_HumanOmmni_forGENESIS.rda'))

## list of individuals in VCF/GDS
write.table(as.data.frame(df[,1]),file='~/scratch/Idstokeep.txt',col.names=FALSE,row.names=FALSE,quote=FALSE)

## Open Rsquare 
txt <- fread('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/rawALL_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.txt',header=F)
txt <- as.data.frame(txt)

names(txt) <- c("CHROM", "POS","SNP","FILTER","AF","MAF","RSQ","STATUS")

# remove R2=
txt$RSQ <- gsub('R2=','', txt$RSQ)
txt$RSQ <- as.numeric(as.character(txt$RSQ))

## NOt all sites were imputed so get the imputes 
imputed <- txt[txt$STATUS=='IMPUTED',]

## get summary of R-square 
print('Summary of R-square')
summary(imputed$RSQ)

## Select sites with with R-square < 0.3
imputed <- imputed[imputed$RSQ < 0.3,c(1,2)]

## NUmber of Site R2<0.3
print('Number of Sites to exclude because R-square < 0.3')
nrow(imputed)

## write list
write.table(imputed,file=paste0(file.path,'SitesExcludeRsq_ALL_Imputed_1_22_FiltDeep_Immune.txt'),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


## Original IMPUTED VCF BEFORE QC
VCF <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/OUTPUT_MIS/VCF/Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.vcf.gz'
NVCF <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/OUTPUT_MIS/VCF/All_QC.Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.vcf.gz'

# list of sites to keep 
ids <- '/users/cvalenci/scratch/Idstokeep.txt'
sites='/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/All_Participants/maf5pct_ALLrawImputed_1_22.txt'
exsites='/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/All_Participants/SitesExcludeRsq_ALL_Imputed_1_22_FiltDeep_Immune.txt'

# Run VCFtools TO CREATE INTERMEDIATE VCF
system(paste0('vcftools',' --gzvcf ',
                      VCF,
                      ' --keep ',ids,
                      ' --maf 0.01 ',
                      ' --exclude-positions ' ,exsites,
                      ' --recode  --recode-INFO-all --stdout ',
                      '| bgzip -c > ',
                      NVCF))

## CREATE GDS FROM NEW VCF (NVCF) INCLUDING DOSAGE
gdsfile <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/All_Participants/All_QC.Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.gds'

seqVCF2GDS(NVCF, gdsfile, verbose=FALSE,fmt.import="DS")


## Reproducibility info
print('Reproducibility information:')
Sys.time()
proc.time()
devtools::session_info()