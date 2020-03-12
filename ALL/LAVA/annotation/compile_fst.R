#!/bin/R

library(ggplot2)
library(data.table)

indir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/'

eur <- system(sprintf('ls %s | grep "weir" | grep "EUR"',indir),intern=TRUE)
afr <- system(sprintf('ls %s | grep "weir" | grep "AFR"',indir),intern=TRUE)

for (pop in c('eur', 'afr')){
	
	if(file.exists(sprintf('%s/%s-merged-fst.tsv',indir,pop))){next}
	first <- TRUE
	for (i in get(pop)){
		dt <- fread(sprintf('%s/%s',indir,i),sep='\t',header=TRUE)
		comparison <- gsub('\\.weir\\.fst','',i)
		dt <- dt[WEIR_AND_COCKERHAM_FST != '-nan']
		colnames(dt) <- c('CHROM','POS',comparison)
		dt <- unique(dt)
		# remove any rows with duplicate positions 
		dt[,snp := paste0(CHROM, '_', POS)]
		remove <- dt[duplicated(snp),snp]
		print(head(remove))
		print(length(remove))
		dt <- dt[!snp %in% remove]
		dt[,snp := NULL]
		if(first){
			master <- dt
			first <- FALSE
		} else {
			master <- merge(master, dt, by=c('CHROM','POS'))
		}
	}
	write.table(master, sprintf('%s/%s-merged-fst.tsv',indir,pop),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
	print(head(master))
}

get_chrom <- function(x){ 
	return(as.integer(gsub('chr','',unname(unlist(strsplit(x,'_')))[1])))
}

get_pos <- function(x){ 
	return(as.integer(unname(unlist(strsplit(x,'_')))[2]))
}

if(!file.exists(paste0(indir, '/snps.RData'))){
	
	snps <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/snps.txt',sep='\t',header=FALSE)
	#snps <- fread(cmd='zcat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.all.snps.txt.gz', header=FALSE, sep='\t')
	colnames(snps) <- 'variant_id'

	snps[,CHROM := sapply(variant_id, get_chrom)]
	snps[,POS := sapply(variant_id, get_pos)]

	save(snps, file=paste0(indir, '/snps.RData'))
} else {
	load(paste0(indir, '/snps.RData'))
}

if(!file.exists(paste0(indir, '/merged-all_no_max.RData'))){
	afr <- fread(paste0(indir,'/afr-merged-fst.tsv'), sep='\t', header=TRUE)
	eur <- fread(paste0(indir,'/eur-merged-fst.tsv'), sep='\t', header=TRUE)
	merged <- merge(afr, eur, by=c('CHROM','POS'), all=TRUE)
	save(merged, file=paste0(indir, '/merged-all_no_max.RData'))
} else {
	load(paste0(indir, '/merged-all_no_max.RData'))
}
merged[ , (colnames(merged)) := lapply(.SD, as.numeric), .SDcols = colnames(merged)]
# get max Fst per row
mycols = colnames(merged)[3:ncol(merged)]
merged[,max_fst := max(abs(.SD),na.rm=TRUE), by=c('CHROM','POS'), .SDcols = mycols]
merged[,which_max_fst := colnames(.SD)[which.max(abs(.SD))], by=c('CHROM','POS'), .SDcols = mycols]

merged <- merged[,.(CHROM, POS, max_fst, which_max_fst)]
# merge with variant_id
master <- merge(merged, snps, by=c('CHROM','POS'))
print(head(master))
save(master, file=paste0(indir, '/master_all.RData'))


# if(!file.exists(paste0(indir, '/master_afr.RData'))){
# 	afr <- fread(paste0(indir,'/afr-merged-fst.tsv'), sep='\t', header=TRUE)
# 	afr[ , (colnames(afr)) := lapply(.SD, as.numeric), .SDcols = colnames(afr)]
# 	# get max Fst per row
# 	mycols = colnames(afr)[3:ncol(afr)]
# 	afr[,max_fst := max(abs(.SD)), by=c('CHROM','POS'), .SDcols = mycols]
# 	afr[,which_max_fst := colnames(.SD)[max.col(.SD)], .SDcols = mycols]

# 	afr <- afr[,.(CHROM, POS, max_fst, which_max_fst)]
# 	# merge with variant_id
# 	master_afr <- merge(afr, snps, by=c('CHROM','POS'))
# 	print(head(master_afr))
# 	save(master_afr, file=paste0(indir, '/master_afr.RData'))
# }
# warnings()

# if(!file.exists(paste0(indir, '/master_eur.RData'))){
# 	eur <- fread(paste0(indir,'/eur-merged-fst.tsv'), sep='\t', header=TRUE)
# 	eur[ , (colnames(eur)) := lapply(.SD, as.numeric), .SDcols = colnames(eur)]
# 	# get max Fst per row
# 	mycols = colnames(eur)[3:ncol(eur)]
# 	eur[,max_fst := max(abs(.SD)), by=c('CHROM','POS'), .SDcols = mycols]
# 	eur[,which_max_fst := colnames(.SD)[max.col(.SD)], .SDcols = mycols]

# 	eur <- eur[,.(CHROM, POS, max_fst, which_max_fst)]
# 	# merge with variant_id
# 	master_eur <- merge(eur, snps, by=c('CHROM','POS'))
# 	print(head(master_eur))
# 	save(master_eur, file=paste0(indir, '/master_eur.RData'))
# } 
# warnings()

# outputs used in compile_coloc.R

