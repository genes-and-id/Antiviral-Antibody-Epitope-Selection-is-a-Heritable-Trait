# SampleManifest_EBV_MANTRA_HLA.R

# TO RUN
# qrsh -l mem_free=10G,h_vmem=10G
# module load R
# Rscript SampleManifest_EBV_MANTRA_HLA.R > SampleManifest_EBV_MANTRA_HLA_log.txt 2>&1


## GENERATES SAMPLE MANIFEST FOR EBV MANTRA ONLY HLA REGION 
## IT CONTAINS INFORMATION FOR ALL GROUPS (VRC-EUR,UKTWINS,VRC-AFR)

# Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(data.table)
data.table::setDTthreads(7)
options(stringsAsFactors=FALSE)

file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/'


## VRC-EUR files 
df <- data.frame(FILES=dir('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA',pattern='mantra.dat',full.names=TRUE))
df$ID <- gsub('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/|_HLA_mantra.dat','',df$FILES) 
df$NAME <- gsub('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/','',df$FILES)
write.table(df,file=paste0(file.path,'SampleManifest_EBV_MANTRA_HLA.txt'),col.names=FALSE,row.names=FALSE,quote=FALSE)


## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()