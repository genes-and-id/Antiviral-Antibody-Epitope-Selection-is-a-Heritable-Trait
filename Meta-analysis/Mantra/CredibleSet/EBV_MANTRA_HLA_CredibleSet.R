# EBV_MANTRA_HLA_CredibleSet.R

## MANTRA ANALYSIS:
## GENERATE 95% CREDIBLE SETS FROM 'EBV_MANTRA_HLA.sh'
## how to run : qsub EBV_MANTRA_HLA_CredibleSet.sh


# Load packages
library(plyr) 
library(dplyr)
library(tidyr)
library(data.table)
data.table::setDTthreads(7) 
options(stringsAsFactors=FALSE)

file.path <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Assoc/2019_07_19/EBV_PHENOS/Imputed/MANTRA/HLA_Results/'

## CREDIBLE SET
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791416/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223059/
# For each SNP in a region, the posterior probability that this SNP is driving the region’s association signal 
# was calculated by dividing the SNP’s BF by the summation of the BFs of all SNPs in the REGION.
# REGION: 500KB +/- top snp
# Top snp: lowest BF
# define region from top-500KB to top+500KB
# pull variants in defined region
# order by BF (decreasing=T)
# calculate posterior for each snp as BF-snp/SUM BF region
# credible set 95%: all ranked variants till posterior Cumsum 0.95

## OPEN RESULTS 
df <- fread(paste0(file.path,'EBV_MANTRA_HLA_Results.txt'),header=T)

names(df)[1] <- 'Peptide'


## List of Peptides 
peptide <- as.list(unique(df$Peptide))
ID <- unique(unlist(peptide))

getCredibleSets <- function(peptide){
	temp <- df[df$Peptide==peptide,]
	temp$BF <- 10^temp$log10Bayes
	topost <- as.numeric(temp[temp$log10Bayes==max(temp$log10Bayes),'pos'])
	min.p <- topost-500000
	max.p <- topost+500000
	temp <- temp[temp$pos>min.p& temp$pos<=max.p,]
	temp <- temp[order(temp$BF,decreasing=T),]
	temp$index <- 1:nrow(temp)
	temp$Posterior <- temp$BF/sum(temp$BF)
	temp$CumPost <- cumsum(temp$Posterior)
	temp <- temp[temp$CumPost<=0.95,]
	temp <- temp[order(temp$pos),-c(1,9)]
}

# Extract results from each GWAS
allresults <- lapply(peptide,getCredibleSets)
names(allresults) <- ID

# make dataframe 
dfall <- ldply(allresults,data.frame)
names(dfall)[1] <- 'Peptide'

## write
write.table(dfall,file=paste0(file.path,'EBV_MANTRA_HLA_CredibleSet.txt'),row.names=FALSE,quote=FALSE)

## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()