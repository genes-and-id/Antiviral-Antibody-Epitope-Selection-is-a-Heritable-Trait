# EBV_Inputs_HLA_mantra.dat.R


## MANTRA ANALYSIS:
## Generate Input files mantra.dat ONLY CHROMOSOME-6 HLA
## HOW TO RUN: qsub EBV_Inputs_HLA_mantra.dat.sh 

# HLA region
#Assembly:
#GRCh37
#Location:
#chr6:28,477,797-33,448,354

# Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(R.utils)
library(data.table)
data.table::setDTthreads(7) # From Leo:to avoid using too many cores
options(stringsAsFactors=FALSE)

# Create directory for ANALYSIS
main_dir <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/'
sub_dir <- 'MANTRA'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}


file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/'
file.data <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/'

## Open main table EBV. This cotains the peptides 
df <- fread(paste0(file.data,'Main.EBV.Table.2020-08-18.csv'))

## GEt the variants from the HLA region
hla <- as.data.frame(df[df$chr==6,])
peptides <- as.list(unique(hla$Peptide))
rm(hla,df)


# Get ID for job 
iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))
## set array job
file <- peptides[[iscen]]
# ID ?
#ID <- unlist(file)

## Extract that data in the following format
# for each variant:
# rsid,chr,pos,effect(alt),other(ref)
# 0/1,n,AF,beta,se

## HLA Region in chr 6 
min <- 28477797
max <- 33448354


## EUR: SELECT SITES BASED ON HLA, FREQ, AND P <0.05
EUR.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/EUR_Genesis/'
EUR <- fread(paste0(EUR.path,file,'.assoc'))
EUR <- EUR[EUR$chr==6,] # select chromosome 
EUR <- EUR[EUR$pos>=min & EUR$pos<=max,] # ONLY HLA REGION
EUR <- EUR[EUR$Score.pval<=0.05,] # ONLY SITES WITH P <= 0.05 
EUR <- as.data.frame(EUR[EUR$freq>=0.05,c('SNP','chr','pos','alt','ref','n.obs','freq','Est','Est.SE')])


## UK
uk.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/Massimo/TwinsUK_QC_2020/'
uk <- fread(paste0(uk.path,'TWINSUK.',file,'.20191112_QCs+RS.txt.result.txt.gz'))
uk <- uk[uk$CHR==6, ]
uk <- uk[uk$pos>=min & uk$pos<=max,]
uk <- as.data.frame(uk[,c('SNP','CHR','BP','effect_allele','other_allele','n','EAF','beta','se')])
names(uk) <- c('SNP','chr','pos','alt','ref','n.obs','freq','Est','Est.SE') # to mach VRC

## AFR 
afr.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/AFR_Genesis/'
afr <- fread(paste0(afr.path,file,'.assoc'))
afr <- afr[afr$chr==6,]
afr <- afr[afr$pos>=min & afr$pos<=max,]
afr <- as.data.frame(afr[ ,c('SNP','chr','pos','alt','ref','n.obs','freq','Est','Est.SE')])

gc()

## Prepapare data 
df <- data.frame( )


for(i in 1:nrow(EUR)){
	snp <- EUR[i,1:5] # get row i, 5 element 
	w <- as.character(EUR[EUR$SNP %in% snp$SNP ,6:9]) # 4 element based on condition presence Y/N
	w <- c('1',w)
	if(nrow(uk[paste0(uk$SNP,uk$ref,uk$alt) %in% paste0(snp$SNP,snp$ref,snp$alt) ,6:9]) >0){
	m <- as.character(uk[paste0(uk$SNP,uk$ref,uk$alt) %in% paste0(snp$SNP,snp$ref,snp$alt) ,6:9]) # 4 element based on condition presence Y/N
	m <- c('1',m)
	} else{
	m <- rep(0,5)
	} # add 1 indicating presence in data
	if(nrow(afr[afr$SNP %in% snp$SNP ,6:9]) >0){
	a <- as.character(afr[afr$SNP %in% snp$SNP ,6:9]) # 4 element based on condition presence Y/N
	a <- c('1',a)
	} else{
	a <- rep(0,5)
	}
	temp <- rbind(as.character(snp[1,]),w,m,a) # EUR,UKtwins,AFR
	df <- rbind(df,temp)
	rm(snp,w,m,a,temp)
	gc()
}


## Write file
## Extract that data in the following format 
# for each variant:
# rsid,chr,pos,effect(alt),other(alt)
# 0/1,n,AF,beta,se
write.table(df,file=paste0(file.path,file,'_HLA','_mantra.dat'),col.names=FALSE,row.names=FALSE,quote=FALSE,sep=" ")

## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()