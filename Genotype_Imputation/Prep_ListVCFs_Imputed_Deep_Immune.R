# Prep_ListVCFs_Imputed_Deep_Immune.R

#qrsh -l mem_free=20G,h_vmem=20G
#module load R

# To run:
## Rscript Prep_ListVCFs_Imputed_Deep_Immune.R > Prep_ListVCFs_Imputed_Deep_Immune_log.txt 2>&1

## Library
library(plyr)
library(devtools)
options(stringsAsFactors=FALSE)

filepath <- '/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/'

# files are located in
files <- data.frame(file=dir('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/OUTPUT_MIS/VCF',pattern='dose.vcf.gz',full.names=TRUE),
	Name=dir('/dcl01/leased/pduggal/UK_TwinData/Deep_Immune_HumanOmni5-4v1-1/CIDR/for_priya_aspera/gwas/Imputation/OUTPUT_MIS/VCF',pattern='dose.vcf.gz',full.names=FALSE))

# remove everything after '.'
files$Index <- gsub("\\..*","",files$Name)
# remove chr
files$Index <- gsub('chr','',files$Index)

# order by Index
files <- files[order(as.numeric(files$Index)), ]

# write.file
write.table(as.data.frame(files[,2]),file=paste0(filepath,'List_vcf.txt'),col.names=FALSE,row.names=FALSE,quote=FALSE)

## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()