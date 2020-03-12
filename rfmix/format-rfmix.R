#!/usr/bin/env Rscript
# Nicole Gay

library(data.table)
options(stringsAsFactors = FALSE)
options(scipen=10000)

# Convert VCFs (admixed individuals and reference panels) to required input format for RFMix

args <- commandArgs(trailingOnly=TRUE)
subjects <- args[1] # text file with admixed subject IDs, one per line, used to filter VCF with admixed subjects
tdictfile <- args[2] # 1000 Genomes metadata file that allows me to select individuals by population
gtexfile <- args[3] # VCF that includes admixed individuals, single chromosome
tgfile <- args[4] # VCF that includes reference panels (1000 Genomes, in this case), single chromosome
outdir <- args[5] # output directory 
chrom <- paste0('chr',args[6]) # specify a single chromosome. run this script in parallel for all chromosomes

########################################################################################
## Import phenotype tables

# GTEx 
# RACE column:
#	1 = Asian
# 	2 = Black or African American
# 	3 = White
# 	4 = American Indian or Alaska Native
# 	99 = Unknown
# 	98 = Not Reported
subjs <- fread(subjects, sep='\t', header=FALSE)
colnames(subjs) <- "ID"
# gpheno <- gpheno[grepl("1|2",gpheno$RACE),]

# 1000 Genomes (Sample, Population)
tdict <- fread(tdictfile, sep='\t', header=TRUE)
# filter for desired populations
# 	CHB Han Chinese in Beijing, China + JPT Japanese in Tokyo, Japan
# 	CEU Utah Residents (CEPH) with Northern and Western European Ancestry
# 	YRI Yoruba in Ibadan, Nigeria + 3 other West African populations
tdict <- tdict[grepl("CHB|JPT|CEU|YRI|GWD|MSL|ESN",Subpop)]

########################################################################################
## filter and format VCFs 

system("echo 'Filtering genotype files...'")
## GTEx FILTERING
# read in admixed VCF
gtex <- fread(cmd=sprintf("zcat %s", gtexfile), sep = '\t', header=TRUE)
# filter subjects by subject IDs in "subjects" input file 
m <- na.omit(match(subjs[,ID],colnames(gtex)))
# keep only columns I want
cols <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',colnames(gtex)[m])
gtex <- gtex[,..cols]
# remove variants that didn't pass filter
gtex <- gtex[FILTER=="PASS"]

## 1000 GENOMES FILTERING
# read in VCF
tg <- fread(cmd=sprintf("zcat %s",tgfile), sep = '\t', header=TRUE)
# order tg columns by class labels; filter subjects
asn <- tdict[Subpop=="CHB"|Subpop=="JPT",Subject] #1
eur <- tdict[Subpop=="CEU",Subject] #2
afr <- tdict[Subpop=="YRI"|Subpop=="GWD"|Subpop=="MSL"|Subpop=="ESN",Subject] #3
o <- c(asn,eur,afr)
m <- na.omit(match(o,colnames(tg)))
# keep only columns I want
cols <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',colnames(tg)[m])
tg <- tg[,..cols]
# filter out multiallelic variants (only applies to 1000 Genomes)
tg <- tg[ALT!=',']
tg <- tg[FILTER=="PASS"]

########################################################################################
## filter down to subset of variants found in both files 

system("echo 'Finding overlap...'")
# select only shared variants 
colsd <- c('ID','QUAL','FILTER','INFO','FORMAT')
tg[,(colsd) := NULL]
gtex[,(colsd) := NULL]
shared <- merge(gtex,tg,by=c("#CHROM","POS","REF","ALT"))
print(paste("N shared variants:", nrow(shared), sep=' '))
if (nrow(shared)==0){
	print(paste0("No shared variants. Aborting execution for ",chrom))
	quit()
}

########################################################################################
## format data to meet RFMix input specs 

system("echo 'Formatting data for RFMix...'")

# save these columns for later
snps <- shared[,list(`#CHROM`,POS,REF,ALT)]
colsd <- c('#CHROM','POS','REF','ALT')
shared[,(colsd) := NULL]

# save vector of subject names (correspond to columns of RMFix input matrix)
subjids <- colnames(shared)
subdf <- data.frame(ID=subjids)
write.table(subdf, file=paste0(outdir, "/subject-ids-",chrom), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)

# add numeric labels (classes) to all of the populations
# 	0 = admixed (all GTEx)
# 	1 = CHB|JPT (Asian)
#	2 = CEU (European)
# 	3 = YRI (African)
get_pop <- function(x){
	if(x %in% tdict[,Subject]){
		if(grepl("CHB|JPT",tdict[Subject==x,Subpop])){
			return(1)
		}else if(grepl("CEU",tdict[Subject==x,Subpop])){
			return(2)
		}else{
			return(3)
		}
	} else {
		return(0)
	}
}
subdf$label <- sapply(subdf$ID, get_pop)

# collapse alleles
system("echo 'Collapsing alleles...'")
write.table(shared, paste0(outdir, "/shared-tmp-", chrom,".txt"), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
system(sprintf('sed -i -e "s/|//g" -e "s/	//g" %s',paste0(outdir, "/shared-tmp-", chrom,".txt")))

write.table(snps, paste0(outdir, "/snps-tmp-", chrom,".txt"), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
# remove rows with missing data
system(sprintf('paste %s %s > %s', paste0(outdir, "/snps-tmp-", chrom,".txt"), paste0(outdir, "/shared-tmp-", chrom,".txt"), paste0(outdir, "/master-tmp-", chrom,".txt")))
system(sprintf('sed -i "/\\./d" %s', paste0(outdir, "/master-tmp-", chrom,".txt")))

# save snp_locations in bed format - THIS IS AN INTERMEDIATE INPUT FOR RFMIX
snp_locations <- fread(paste0(outdir, "/master-tmp-", chrom,".txt"), sep='\t', header=FALSE, select=c(1,2), col.names=c("#CHROM","POS"))
snp_locations[,START:=POS-1]
colnames(snp_locations) <- c("CHR","START","END")
write.table(snp_locations, file=paste0(outdir, "/snp-bed-", chrom), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)

# save the binary matrix - THIS IS AN INPUT FOR RFMIX
system("echo 'Writing out binary alleles'")
system(sprintf('cut -f 5 %s > %s', paste0(outdir, "/master-tmp-", chrom,".txt"), paste0(outdir,"/alleles-", chrom)))

# save class vector for RFMix input - THIS IS AN INPUT FOR RFMIX
classvector <- rep(subdf$label, each = 2) # needs to be duplicated: 2 haplotypes for each subject
classes <- t(data.frame(classvector))
system("echo 'Saving subject class file'")
write.table(classes, file=paste0(outdir, "/classes-", chrom), col.names=FALSE, row.names=FALSE, sep=' ', quote=FALSE)

# remove tmp files
system(sprintf('rm %s %s %s', paste0(outdir, "/snps-tmp-", chrom,".txt"), paste0(outdir, "/shared-tmp-", chrom,".txt"), paste0(outdir, "/master-tmp-", chrom,".txt")))
