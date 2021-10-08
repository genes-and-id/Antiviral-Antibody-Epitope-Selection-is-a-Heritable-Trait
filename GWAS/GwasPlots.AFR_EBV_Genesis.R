# GwasPlots.AFR_EBV_Genesis.R


## TO RUN: qsub GwasPlots.AFR_EBV_Genesis.sh

## RESULTS FROM GWAS.EBV.AFR_Genesis_2019-07-19.R
# this code take the GWAS results and produced Manhattan plots and Q-Q plots
# for each independant GWAS

## libraries
library(qqman)
library(data.table)
data.table::setDTthreads(7)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo') # solve problem with X11
Sys.setenv("DISPLAY"=":0")


# Create directory for ANALYSIS
main_dir <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis/'
sub_dir <- 'Plots'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

## Define IDs
# list data frames 
all_df <- list.files(path="/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis",pattern=".assoc",full.names=TRUE)
all_df <- as.list(all_df)

## Get ID for job 
iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## set array job
file <- all_df[[iscen]]

## Get ID for df 
ID <- gsub('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis/|.assoc','',file)

## data
# open data 
df <- fread(file)
df <- as.data.frame(df)
df <- df[df$freq>=0.05,] ## REMOVE SITES WITH freq < 0.05

# take snp,chr,pos,bs,p-value
temp <- df[,c(2,1,3,11)]
# remove NA's
temp <- temp[!is.na(temp$Score.pval), ]

## Inflation factor lambda 

# get p values
pvalues <- temp$Score.pval
chisq <- qchisq(1-pvalues,1)
# calculate lambda & round value
lambda <- round(median(chisq)/qchisq(0.5,1),3)

# setwd for the plots
setwd('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis/Plots')

## Manhattan plot
png(filename=paste0('Manhattan.EBV_AFR_',ID,'.png'),width=1800,height=800,units="px",pointsize=12,bg="white")
manhattan(temp, chr = "chr", bp = "pos", p = "Score.pval",
			snp = "variant.id", main=paste('EBV-AFR',ID), 
			ylim=c(0,10),chrlabs = NULL,suggestiveline = -log10(1e-05),
			genomewideline = -log10(5e-08),highlight = NULL,logp = TRUE,
			cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)		
dev.off()

## Q-Q plot
png(filename=paste0('qqplot_EBV_AFR_',ID,'.png'),width=1800,height=800,units="px",pointsize=12,bg="white" )
qq(temp$Score.pval,main=paste('EBV-AFR',ID,'Lambda=',lambda),
cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5 )
dev.off()


## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()