# Author: Nicole Gay
# Winter 2018
# Inputs: Concatenated rows for pca, platform, and sex covariates from all GTEx v8 covariate files; PEER covariates
# Output: all_covariates file
# Called from covariate-format.sh

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
peer <- args[2]
outfile <- args[3]

newcov <- read.delim(infile, sep='\t', header=FALSE, stringsAsFactors=FALSE)
remove <- duplicated(as.character(newcov[1,]))
head(newcov)
newcov[,remove] = NULL
rownames(newcov) <- c("ID","pcr","platform","sex")

cov <- read.delim(peer, header=FALSE, stringsAsFactors=FALSE)

d <- newcov[,match(as.character(cov[1,]),as.character(newcov[1,]))]
d <- d[-1,]
colnames(cov) <- c(1:ncol(cov))
colnames(d) <- c(1:ncol(cov))
new <- rbind(cov,d)

write.table(new, outfile, sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
