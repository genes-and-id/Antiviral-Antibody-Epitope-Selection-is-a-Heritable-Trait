# GWAS.EBV.EUR_Genesis_2019-07-19.R

## ARRAY JOB
# TO RUN : qsub GWAS.EBV.EUR_Genesis_2019-07-19.sh 


## RUNS GWAS ON EBV PEPTIDES 
## IT USES IMPUTED GENOTYPES 
## INPUTED GENOTYPES: final_EA_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.gds
## OUTPUT: Pepetide.assoc FOR EACH PEPTIDE

# Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(SeqVarTools)
library(SeqArray)
library(Rsamtools)
library(data.table)
data.table::setDTthreads(7) # From Leo:to avoid using too many cores
options(stringsAsFactors=FALSE)


# Create directory for ANALYSIS
main_dir <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/'
sub_dir <- 'EUR_Genesis'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/'
file.data <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/All_Participants/'


## LOAD PHENO: THIS HAS ONLY FIRST 10 PCS
load(paste0(file.data,'EBV_pheno2019-07-19_HumanOmmni_forGENESIS.rda'))
df <- df
# 631 
## first column: scanID
# PCs from 2:11,  PC1,PC2..
# Pheno from 12:108,  Peptide

## Make a list of Pepetides 
all_peptides <- as.list(names(df)[12:ncol(df)])
#all_peptides <- as.list(names(df)[12:13]) # testing

## Get ID for job 
iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## set array job for each PEPTIDE
peptide <- all_peptides[[iscen]]

## MAKE PHENO FOR GWAS: scanID,PCs,then peptide
#vars <- names(df)[1:11] # 1-10 PCs but if not PCs then use line below
vars <- names(df)[1]
vars <- c(vars,peptide)
phenogwas <- df[,colnames(df) %in% vars]
names(phenogwas)[1] <- 'sample.id'
# recode peptide 
# it was : case=2 to 1 ; control=1 to 0; NA=0 to 'NA'
pp <- names(phenogwas)[ncol(phenogwas)]
phenogwas[pp] <- lapply(phenogwas[pp],function(x) ifelse(x==2,1,ifelse(x==1,0,NA))  )
## ID for PEPTIDE
ID <- unlist(names(phenogwas)[ncol(phenogwas)])
## fix name for GENESIS (silly)
names(phenogwas)[ncol(phenogwas)] <- 'pheno'

## Select Individuals
# this ok, it will introduce NA's but it will match gds 
individuals <- read.table('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/Plink_Files/final_EA_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.fam',header=F)
addids <- data.frame(sample.id=individuals[,1])

phenogwas <- merge(phenogwas,addids,by='sample.id',all.y=TRUE)

# Checking
print('Peptide name')
vars[length(vars)]
ID


## READING GDS FILE
geno <- seqOpen('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/Genesis_EUR/final_EA_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.gds')

# make SeqVarIteratetor object
seqData <- SeqVarData(geno, sampleData=AnnotatedDataFrame(phenogwas))

# Define block SNPs using default variantBlock=10000
iterator <- SeqVarBlockIterator(seqData)


## GWAS: Generalized linear mix model (GLMM)

# define covariates (PCs) : NO FOR EUR

## define null model with 
# , covars = cov
# cov.mat = myGRM,
# ,sample.id=individuals

nullmod <- fitNullModel(iterator, outcome = 'pheno' , family = 'binomial')

## RUN ASSOCIATION USING DOSAGE: this return a data frame
assoc <- assocTestSingle(iterator, null.model = nullmod,imputed=TRUE)



## Add SNPID 
id <- fread('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Clean.Data/FINAL_QC/snpID.final_EA_Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.vcf.txt')
id <- as.data.frame(id)
names(id) <- c('chr','pos','SNP','ref','alt')
# add index for variant.id
id$variant.id <- 1:nrow(id)
# merge using chr,pos,variant.id 
assoc <- merge(assoc,id,by=c('chr','pos','variant.id'))
# re-order columns
assoc <- assoc[,c(1,2,15,16,17,5:14)]
## sort 
assoc <- assoc[order(assoc$chr,assoc$pos),]

## write assoc
write.table(assoc,file=paste0(file.path,ID,'.assoc'),row.names=FALSE,quote=FALSE )

## clean 
rm(list=ls())


## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()