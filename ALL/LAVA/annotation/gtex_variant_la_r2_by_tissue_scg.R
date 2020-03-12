# 9/23/2019

# find VE between GTEx v8 lead SNP genotype and LA, accounting for individuals in each tissue
# built off of r2_gtex_lead_snps.R

library(data.table)
library(ggplot2)
library(OneR)
library(plyr)

args <- commandArgs(trailingOnly=TRUE)
tissue <- args[1] # do this for all tissues in GTEx

outdir='/labs/smontgom/nicolerg/LAVA/r2_geno_la/results/'
indir='/labs/smontgom/nicolerg/LAVA/r2_geno_la/data/'

load(sprintf('%s/%s-data-for-regressions.RData',indir,tissue)) # test_snps, vcf_filt, sub_local

# create a progress bar
test_snps <- unique(test_snps)
total <- length(test_snps)
count = 0
pb <- txtProgressBar(min = 0, max = total, style = 3)

outfile <- file(paste0(outdir,'/',tissue,'_leadsnp.r2.txt'), 'w')
cat("variant_id\tlocal_rsq\n", file=outfile)

try_lm <- function(f, df) {
    out <- tryCatch(
        {
			summary(lm(f, data=df))$r.squared
        },
        error=function(cond) {
            return(NA)
        }
    )    
    return(out)
}

for (v in test_snps){
	
	# skip SNPs not in this file
	if (!v %in% colnames(vcf_filt)){
		next
	}

	g <- vcf_filt[,c('SUBJID', v), with=FALSE]
	l <- sub_local[SNP_ID == v]
	m <- merge(g, l, by='SUBJID')
	f <- as.formula(sprintf('%s ~ AFR + ASN',v))

	rsq <- try_lm(f, m)

	string <- paste(v, rsq, sep='\t')
	cat(paste0(string,'\n'), file=outfile, append=TRUE)

	count <- count + 1
	setTxtProgressBar(pb, count)
}
close(pb)
close(outfile)
