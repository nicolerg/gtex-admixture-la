library(SNPRelate)
library(data.table)
library(ggplot2)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)
pop <- args[1]
chr <- args[2]

peerdir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/peer/"
covdir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/cov/"

if (pop=="EUR"){
	geuv <- read.delim(paste0(covdir,"geuvadis-eur-cov.tsv"), sep='\t', header=TRUE, stringsAsFactors=FALSE)
	peer <- read.delim(paste0(peerdir,"GV-EUR-chr",chr,".PEER_covariates.txt"), sep='\t', header=TRUE, stringsAsFactors=FALSE)

} else if (pop=="AFR"){
	geuv <- read.delim(paste0(covdir,"geuvadis-afr-cov.tsv"), sep='\t', header=TRUE, stringsAsFactors=FALSE)
	peer <- read.delim(paste0(peerdir,"GV-AFR-chr",chr,".PEER_covariates.txt"), sep='\t', header=TRUE, stringsAsFactors=FALSE)
}

geuv <- geuv[,c("Sample","Gender")] 

# merge Geuvadis covariates and PEER factors
ids <- peer$ID 
peer$ID <- NULL
subjs <- colnames(peer)
peer <- data.frame(t(peer),stringsAsFactors=FALSE)
colnames(peer) <- ids
peer$SUBJECT <- subjs
merged <- merge(geuv, peer, by.x="Sample", by.y="SUBJECT")

# format merged covariates table: covariate ID x subject
subjs <- merged$Sample
merged$Sample <- NULL
cols <- colnames(merged)
merged <- data.frame(t(merged),stringsAsFactors=FALSE)
colnames(merged) <- subjs
merged$ID <- cols
merged <- merged[,c("ID",subjs)]
rownames(merged) <- NULL
write.table(merged, paste0(covdir,pop,"-chr",chr,".covariates.tsv"), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

