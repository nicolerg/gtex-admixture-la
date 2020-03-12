#!/usr/bin/env Rscript
# Nicole Gay

library(data.table)
options(stringsAsFactors = FALSE)
options(scipen=10000)

# Convert VCFs (admixed individuals and reference panels) to required input format for RFMix

args <- commandArgs(trailingOnly=TRUE)
vcf <- args[1] # merged VCF
tg_meta <- args[2]
chrom <- args[3]
outdir <- args[4] # output directory 

# format inputs for RFMix

vcf <- fread(cmd=sprintf('zcat %s',vcf), sep='\t', header=T, skip="#CHROM")
tg_meta <- fread(tg_meta, sep='\t', header=T)

subjects <- colnames(vcf)[10:ncol(vcf)]

# add numeric labels (classes) to all of the populations
# 	0 = admixed (all GTEx)
# 	1 = CHB|JPT (Asian)
#	2 = CEU (European)
# 	3 = YRI (African)
get_pop <- function(x){
	if(x %in% tg_meta[,Subject]){
		if(grepl("CHB|JPT",tg_meta[Subject==x,Subpop])){
			return(1)
		}else if(grepl("CEU",tg_meta[Subject==x,Subpop])){
			return(2)
		}else{
			return(3)
		}
	} else {
		return(0)
	}
}
classes <- unname(sapply(subjects, get_pop))
classes <- rep(classes, each=2)

# write out "classes" file - THIS IS AN INPUT FOR RFMIX
write(classes, file=sprintf('%s/rfmix_inputs/classes-chr%s',outdir,chrom), ncolumns=length(classes), sep=' ')

# write out "snp_positions" file - THIS IS AN INTERMEDIATE FILE
snps <- vcf[,c(1:3),with=F]
snps[,START := POS - 1]
# #CHROM POS ID START
colnames(snps) <- c('CHROM','STOP','VARIANT_ID','START')
setcolorder(snps, c('CHROM','START','STOP','VARIANT_ID'))
write.table(snps, file=sprintf('%s/rfmix_inputs/snp_positions-chr%s.bed',outdir,chrom), col.names=F, row.names=F, quote=F, sep='\t')

# write out "alleles" file - THIS IS AN INPUT FOR RFMIX
alleles <- vcf[,c(1:9) := NULL]
write.table(alleles, file=sprintf('%s/rfmix_inputs/alleles-chr%s',outdir,chrom), col.names=F, row.names=F, quote=F, sep='')
system(sprintf('sed -i "s/|//g" %s/rfmix_inputs/alleles-chr%s',outdir, chrom))
