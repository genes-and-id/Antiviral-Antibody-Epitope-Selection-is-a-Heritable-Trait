# Inputs_EBV_Genesis_to_Meta.R

## TO RUN: qsub Inputs_EBV_Genesis_to_Meta.sh

## PREPARE INPUT FILES FOR META-ANALYSIS WITH META FOR EBV VRC EURS

## libraries
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(R.utils) # to open gz files
library(data.table)
data.table::setDTthreads(7)
options(stringsAsFactors=FALSE)

# Create directory for ANALYSIS
main_dir <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/'
sub_dir <- 'META'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/META/'

# FILES GENERATED PREVIOUSLY FOR 'METAL'
all_df <- list.files(path='/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/For_Metal',pattern='.gz',full.names=TRUE)

## Match peptides from VRC-EURs to TWINSUK 
Peptides <- data.frame(FILES=dir('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/Massimo/TwinsUK_QC_2020', pattern='txt.result.txt.gz', full.names=FALSE))
Peptides$Name <- gsub('TWINSUK.','',Peptides$FILES)
Peptides$Name <- gsub("\\..*","",Peptides$Name)
Peptides$Name <- paste0('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/For_Metal/forMetal_',Peptides$Name,'.assoc.gz')

## Filter all_df to keep only the common peptides 
all_df <- all_df[all_df %in% Peptides$Name]

# IF missing peptides then remove comment (#) and do: 
#missing_peptides <- data.frame(FILES=dir(file.path,pattern='.assoc.txt.gz'))
#missing_peptides$FILES <- gsub('\\..*','',missing_peptides$FILES)
#missing_peptides$FILES <- paste0('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/For_Metal/forMetal_',missing_peptides$FILES,'.assoc.gz')

#all_df <- all_df[!all_df %in% missing_peptides$FILES]

## Get ID for job 
iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))
## set array job
file <- all_df[[iscen]]

## Get ID for the failed file: 
IDGZ <- gsub('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/For_Metal/forMetal_|.assoc.gz','',file)


## open data 
df <- fread(file)

## Open file with R2 : RSQ we will called info 
anot <- fread('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/All_QC.Sort_Imputed_1_22_FiltDeep_Immune_ibdmafhwe.txt',header=F)
names(anot) <- c("chr", "pos","SNP","FILTER","AF","MAF","info","STATUS")
# remove R2=
anot$info <- gsub('R2=','', anot$info)
anot$info <- as.numeric(as.character(anot$info)) # introduced NA on sites label as 'TYPED_ONLY'
anot <- anot[,c(1:3,7)]

## merge to add R2
temp <- merge(df,anot,by=c('chr','pos','SNP'),all.x=TRUE)
temp$info[is.na(temp$info)] <- 'TYPED_ONLY' # FIX NAs

# fix names for meta:
# SNP=rsid
# Score.pval=P_value
# Est=beta
# Est.SE=se
# allele_A=non-effect allele; IN OUR CASE REF ALLELE (0)
# allele_B=effect allele; ALT ALLELE (1)
#names(temp)[c(3,11,12,13,15,16)] <- c('rsid','P_value','beta','se','allele_A','allele_B')

## META REQUIRES INPUT IN THIS EXACT FORMAT 
## BOTH ORDER OF COLUMNS AND COLUMNS NAMES AND NO EXTRA COLUMNS
temp <- temp[,c(1,3,2,15,16,17,11,12,13)]
# name columns following META FORMAT
names(temp) <- c('chr','rsid','pos','allele_A','allele_B','info','P_value','beta','se')

## write space delimited gzipped txt file 
write.table(temp,file=gzfile(paste0(file.path,IDGZ,'.assoc.txt.gz' )), row.names=FALSE,quote=FALSE,sep=" " )


rm(temp)
gc()

## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()