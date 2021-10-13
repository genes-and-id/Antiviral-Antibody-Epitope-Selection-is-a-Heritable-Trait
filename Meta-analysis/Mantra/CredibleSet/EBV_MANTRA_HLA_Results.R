# EBV_MANTRA_HLA_Results.R


## MANTRA ANALYSIS:
## CREATE FILE WITH ALL RESULTS FROM 'EBV_MANTRA_HLA.sh'
## HOW TO RUN: qsub EBV_MANTRA_HLA_Results.sh


# Load packages
library(plyr) 
library(dplyr)
library(tidyr)
library(data.table)
data.table::setDTthreads(7) 
options(stringsAsFactors=FALSE)

file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/HLA_Results/'


## mantra.out FILES 
manifest <- read.table('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/SampleManifest_EBV_MANTRA_HLA.txt',header=F)
peptides <- manifest$V2
peptides <- as.list(peptides)


## FUNCTION
getresults <- function(filename){
	df <- fread(paste0(file.path,filename,'/mantra.out'), header=F)
	names(df) <- c('SNP','chr','pos','alt','ref','N.Studies','log10Bayes','posteriorP','n','direction')
	df <- df[df$log10Bayes>0,] # filter sites by BF
	temp <- fread(paste0(file.path,filename,'/mantra.bet.out'),header=F,fill=TRUE)
	temp <- temp[,c(1:3,6,7,9,10,12,13)]
	names(temp) <- c('SNP','chr','pos','VRCEUR.Est','VRCEUR.Est.SD','UKT.Est','UKT.Est.SD','VRCAFR.Est','VRCAFR.Est.SD')
	df <- merge(df,temp,by=c('SNP','chr','pos'))
}


## Get ID for each peptide 
ID <- unique(unlist(peptides))


# Extract results from each GWAS
allresults <- lapply(peptides,getresults)
names(allresults) <- ID

# make dataframe 
dfall <- ldply(allresults,data.frame)

## Write file
write.table(dfall,file=paste0(file.path,'EBV_MANTRA_HLA_Results.txt'),row.names=FALSE,quote=FALSE )


## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()