# EBV_White_Genesis_SuppTable.R


## RUN THE BASH SCRIPT: qsub EBV_White_Genesis_SuppTable.sh


# Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(qqman)
library(R.utils) # to open gz files
library(data.table)
data.table::setDTthreads(7) # From Leo:to avoid using too many cores
options(stringsAsFactors=FALSE)


## DEFINE FILE.PATH
file.path <-'/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/White_Genesis/'

# function to extract results.
# selecting sites 
getresults <- function(filename) {
		df <- fread(filename,header=TRUE)
		df <- df[df$freq>=0.05,] ## REMOVE SITES WITH freq.White < 0.05
		df <- df[df$Score.pval<=5e-6, c(1:7,12:14) ] # Extract columns :chr pos SNP ref alt n.obs freq Score.pval Est Est.SE
}

## RESULTS
# Make a list of the plink results
all_df <- list.files(path="/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/White_Genesis",pattern=".assoc",full.names=TRUE)
all_df <- as.list(all_df)
#all_df <- all_df[1:10] # testing 

## Get ID for each peptide 
ID <- unique(unlist(all_df))
ID <- gsub('\\/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/White_Genesis/|.assoc','',ID)


# Extract results from each GWAS
allresults <- lapply(all_df,getresults)
names(allresults) <- ID

# make dataframe 
dfall <- ldply(allresults,data.frame)
# correct name peptide re-order columns 
colnames(dfall)[1] <- 'Peptide'
dfall <- dfall[,c(2:ncol(dfall),1)]

## write this file as csv: THIS IS FOR RAJA
write.csv(dfall,file=paste0(file.path,'EBV_White_Genesis_5e6.csv'),row.names=FALSE )

## write this file for supplemental material 
write.csv(dfall[dfall$Score.pval<=5e-8, ],file=paste0(file.path,'EBV_White_Genesis_5e8.csv'),row.names=FALSE )


## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()