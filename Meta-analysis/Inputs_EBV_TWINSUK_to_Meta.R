# Inputs_EBV_TWINSUK_to_Meta.R

## TO RUN: qsub Inputs_EBV_TWINSUK_to_Meta.sh

## PREPARE INPUT FILES FOR META-ANALYSIS WITH META FOR EBV TWINSUK

## libraries
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(R.utils) # to open gz files
library(data.table)
data.table::setDTthreads(7)
options(stringsAsFactors=FALSE)

# Create directory for ANALYSIS
main_dir <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/Massimo/TwinsUK_QC_2020/'
sub_dir <- 'META'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/Massimo/TwinsUK_QC_2020/META/'

# FILES GENERATED PREVIOUSLY FOR 'METAL'
all_df <- list.files(path='/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/Massimo/TwinsUK_QC_2020',pattern='txt.result.txt.gz',full.names=TRUE)

## Get ID for job 
iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))
## set array job
file <- all_df[[iscen]]

## Get ID for the failed file: 
IDGZ <- gsub('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/Massimo/TwinsUK_QC_2020/TWINSUK.','',file)
IDGZ <- gsub('\\..*','',IDGZ) 

## Open file 
df <- fread(file)

## Rename columns for META
# this file has SNP,rsid,snpid
# we will nename rsid as: original_rsid
# and SNP as rsid 
# other_allele= allele_A: non-effect allele; IN OUR CASE REF ALLELE (0)
# effect_allele=allele_B: effect allele; ALT ALLELE (1)
# p=P_value

names(df)[c(2,4,8,9,14)] <- c('rsid','original_rsid','allele_A','allele_B','P_value')

## write space delimited gzipped txt file 
write.table(df,file=gzfile(paste0(file.path,IDGZ,'.TWINSUK.txt.gz' )), row.names=FALSE,quote=FALSE,sep=" " )

## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()
