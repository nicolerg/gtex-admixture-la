# Author: Nicole Gay
# 17 May 2019
# Inputs: Concatenated rows for genotype, pca, platform, and sex covariates from all GTEx v8 covariate files; PEER covariates
# Output: all_covariates file

library(data.table)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
peer <- args[2]
outfile <- args[3]

cov <- fread(infile, sep='\t', header=TRUE)
peer <- fread(peer, header=TRUE, sep='\t')

cov <- cov[,colnames(cov)[colnames(cov) %in% colnames(peer)], with=FALSE]
new <- rbind(cov, peer)

write.table(new, outfile, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
