# Meta_Manhattan_EUR_Genesis_EBV_Imputed.R

## PLOTS FOR Whites-EBV GWAS USING GENESIS

## ssh -X cvalenci@jhpce01.jhsph.edu -o ForwardX11Timeout=24h
## qrsh -l mem_free=10G,h_vmem=10G

## module load R

## RUN THE BASH SCRIPT: qsub Meta_Manhattan_EUR_Genesis_EBV_Imputed.sh

# Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(qqman)
library(R.utils) # to open gz files
library(data.table)
data.table::setDTthreads(7) # From Leo:to avoid using too many cores
options(stringsAsFactors=FALSE)

options(bitmapType='cairo') 
Sys.setenv("DISPLAY"=":0")

## DEFINE FILE.PATH
file.path <-'/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/'

## FROM Sidak-Nyholt
# For paper and based on h2 (TWINSUK) and prevalence VRC
# 57 GWAS
# 1.036483e-09
# K.indep: 48.24008

## Single trait Genomewide 5e-08

# function to extract results
# selecting sites 
getresults <- function(filename) {
		df <- fread(filename,header=TRUE)
		df <- as.data.frame(df)	
		df <- df[df$freq>=0.05,] ## REMOVE SITES WITH freq.EUR < 0.05
		df <- df[df$Score.pval <=0.05,c('SNP','chr','pos','Score.pval')] # snp,chr,bp,p ALL WE NEED FOR MANHATTAN
}

## RESULTS
# Make a list of the plink results
all_df <- list.files(path="/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis",pattern=".assoc",full.names=TRUE)
all_df <- as.list(all_df)

## Get ID for each peptide 
ID <- unique(unlist(all_df))
ID <- gsub('\\/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/|.assoc','',ID)

# Extract results from each GWAS
allresults <- lapply(all_df,getresults)
names(allresults) <- ID

# make dataframe 
dfall <- ldply(allresults,data.frame)
# correct name peptide re-order columns 
colnames(dfall)[1] <- 'Peptide'
dfall <- dfall[,c(2:ncol(dfall),1)]

## write this file 
write.table(dfall,file=paste0(file.path,'CombinedManhattan.cv.txt'),row.names=FALSE,quote=FALSE)

## write file with only GWAS hits
write.table(dfall[dfall$Score.pval<=5e-08,],file=paste0(file.path,'GWAS_summary_report_CV.txt'),row.names=FALSE,quote=FALSE )

 

## FILTER PEPTIDES BASED ON OVERLAPP BETWEEN UKTWINS AND VRC
overlap <- fread('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/EBV_52_2019-07-19_HeritableTwins2021-02-23.txt',header=T)

## filter
dfall <- dfall[dfall$Peptide %in% overlap$Peptide,]

# No snpsOfInterest
multihits <- as.data.frame(table(dfall[dfall$Score.pval<=1.036483e-09,3 ])) # 
snpsOfInterest <- as.character(unique(multihits[multihits$Freq>1,1])) # 

# Number of Peptides with sig associations
Peptides <- as.data.frame(table(dfall[dfall$Score.pval<=1.036483e-09,5 ])) # 

## ID 
ID <- 'EUR-EBV-Genesis'
n <- as.character(length(unique(dfall$Peptide)))

## Generate Manhattan plot
#,highlight = snpsOfInterest
png(filename=paste0(file.path,'Plots/','Meta-Manhattan.',ID,'.Imputed.png'),width=1800,height=800,units="px",pointsize=14,bg="white")
manhattan(dfall, chr = "chr", bp = "pos", p = "Score.pval",
          snp = "SNP", main=paste(ID,n,'peptides'),
          ymax=10,chrlabs = NULL,suggestiveline = -log10(5e-08),
          genomewideline = -log10(1.036483e-09),logp = TRUE,
          cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
dev.off()

## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()