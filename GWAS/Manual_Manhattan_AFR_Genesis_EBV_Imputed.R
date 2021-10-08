# Manual_Manhattan_AFR_Genesis_EBV_Imputed.R

## PLOTS FOR AFR-EBV GWAS USING GENESIS

## ssh -X cvalenci@jhpce01.jhsph.edu -o ForwardX11Timeout=24h
## qrsh -l mem_free=10G,h_vmem=10G

## module load R

## RUN THE BASH SCRIPT: qsub Manual_Manhattan_AFR_Genesis_EBV_Imputed.sh

# Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(qqman)
library(R.utils) # to open gz files
library(data.table)
data.table::setDTthreads(7) # From Leo:to avoid using too many cores
options(stringsAsFactors=FALSE)

options(bitmapType='cairo') # solve problem with X11
Sys.setenv("DISPLAY"=":0")

## DEFINE FILE.PATH
file.path <-'/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis/'

## FIRST EXTRACT RESULTS FOR META-MANHATTAN
# function to extract results
getresults <- function(filename) {
		df <- fread(filename,header=TRUE)
		df <- as.data.frame(df)	
		df <- df[df$Score.pval <=0.05,c(2,1,3,11)] # snp,chr,bp,p ALL WE NEED FOR MANHATTAN
}

## RESULTS
# Make a list of the plink results
all_df <- list.files(path="/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis",pattern=".assoc",full.names=TRUE)
all_df <- as.list(all_df)

## Get ID for each peptide 
ID <- unique(unlist(all_df))
ID <- gsub('\\/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis/|.assoc','',ID)


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

## THIS IS SIMPLER THAN BELOW 
multihits <- as.data.frame(table(dfall[dfall$Score.pval<=5e-08,3 ]))
snpsOfInterest <- as.character(unique(multihits[multihits$Freq>1,1]))


## ID 
ID <- 'AFR-EBV-Genesis'
n <- as.character(length(unique(dfall$Peptide)))

## Generate Manhattan plot
png(filename=paste0(file.path,'Plots/','Meta-Manhattan.',ID,'.Imputed.png'),width=1800,height=800,units="px",pointsize=14,bg="white")
manhattan(dfall, chr = "chr", bp = "pos", p = "Score.pval",
          snp = "SNP", main=paste(ID,n,'peptides'),
          ymax=10,chrlabs = NULL,suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-08),highlight = snpsOfInterest,logp = TRUE,
          cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
dev.off()

## Reproducibility info
print('Reproducibility Information:')
devtools::session_info()