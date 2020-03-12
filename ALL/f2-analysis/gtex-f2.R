# extract f2 variants in GTEx admixed individuals

library(data.table)
library(tidyr)
options(scipen=999)

args <- commandArgs(trailingOnly=TRUE)
chr <- as.numeric(args[1])

outdir <- "/mnt/lab_data/montgomery/nicolerg/f2-variants/1_gtex-admix-all-f2var"

# master <- fread('/mnt/lab_data/montgomery/nicolerg/rfmix-outputs/afr.ancestry.bed', sep='\t', header=FALSE)
# colnames(master) <- c("chr","start","end","afr","haplotype")
# master[,afr:=NULL]

gtexids <- read.table('/mnt/lab_data/montgomery/nicolerg/rfmix-inputs/gtex-ids', sep='\t', header=FALSE)
ids <- as.character(unname(unlist(gtexids[,1])))

getchr <- function(x){
	return(unname(unlist(strsplit(x, "_"))[1]))
}

vcf <- paste0('/mnt/lab_data/montgomery/nicolerg/filtered/gtex-phased-unfiltered-hg19/gtex.chr',chr,'.vcf.gz')
geno <- fread(input=sprintf("zcat %s", vcf), sep='\t', header=TRUE)

# filter subjects
var <- geno[,10:ncol(geno)] # genotypes only
del <- setdiff(colnames(var),ids)
var[,(del):=NULL]

# replace ".|." with "3|3"
var <- data.frame(var, check.names=FALSE, stringsAsFactors=FALSE)
var <- data.frame(apply(var, 2, function(x) gsub(".|.","3|3",x,fixed=TRUE)),check.names=FALSE, stringsAsFactors=FALSE)

# convert individuals to haplotype
for (col in colnames(var)){
	var <- separate(var, get(col), into=c(paste0(col,'_A'),paste0(col,'_B')))
}
# convert to numeric 
var <- data.frame(apply(var, 2, function(x) as.numeric(x)), check.names=FALSE, stringsAsFactors=FALSE)
# replace 3 with NA
var[var == 3] <- NA

# get row sums; add back variant info; only keep rows with a row sum of 2
keep <- rowSums(var, na.rm=TRUE) == 2
var$ID <- geno[,ID]
var$POS <- geno[,POS]
var$REF <- geno[,REF]
var$ALT <- geno[,ALT]
var$CHROM <- sapply(var$ID, getchr)
var <- var[keep,]

f2 <- copy(var)
f2 <- data.table(f2[,c("ID","POS","REF","ALT","CHROM")])

subj1 <- rep(NA, nrow(var))
subj2 <- rep(NA, nrow(var))
for(i in 1:nrow(var)){
	vec <- unname(unlist(var[i,]))
	subjs <- colnames(var)[which(vec == 1)]
	if(length(subjs)==2){
		subj1[i] <- subjs[1]
		subj2[i] <- subjs[2]
	} else {
		print("Error: Not an f2 variant",quote=FALSE)
	}
}

f2[,SUBJ1:=subj1]
f2[,SUBJ2:=subj2]
colnames(f2)[colnames(f2)=="POS"] <- "START"
f2[,END:=f2[,START]+1]

# remove variants that are not on this chromosome (don't know why this is happening)
f2[ CHROM == paste0("chr",chr) ]

setcolorder(f2, c("CHROM", "START", "END", "ID", "REF", "ALT", "SUBJ1", "SUBJ2"))
write.table(f2, file=paste0(outdir,"/f2variants-chr",chr,".bed"), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
