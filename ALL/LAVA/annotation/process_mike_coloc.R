## DEPRECATED: for FINEMAP

# get high-confidence colocalizations specifically in genes that are eGenes in both methods (1e-5) with different lead SNPs

library(data.table)
library(ggplot2)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

# match allpairs files with filtered coloc results 

# do it separately for each tissue 

coloc <- fread('/users/mgloud/projects/nicoles-local-ancestry/output/all_finemap_results.txt', sep='\t', header=TRUE)
coloc[,`-log_eqtl_pval`:=NULL]
coloc <- coloc[`-log_gwas_pval` > 4]

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

get_tissue <- function(x){
	return(unname(unlist(strsplit(x, '_filtered')))[1])
}

get_ref_snp <- function(x){
	y <- paste(unlist(unname(strsplit(x, '_')))[1:2], collapse='_')
	y <- gsub("chr","",y)
	return(y)
}

get_method <- function(x){
	if(grepl('local',x)){
		return('LAVA')
	}
	if(grepl('global',x)){
		return('global')
	}
}

coloc[,tissue := sapply(eqtl_file, get_tissue)]
# feature == gene_id

for (t in tissues){
	
	c <- coloc[tissue == t]

	global_filt <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/allpairs_difflead/',t,'.filtered.global.allpairs.tsv.gz')
	local_filt <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/allpairs_difflead/',t,'.filtered.local.allpairs.tsv.gz')
	gfilt <- fread(cmd=sprintf('zcat %s', global_filt), sep='\t', header=TRUE)
	lfilt <- fread(cmd=sprintf('zcat %s', local_filt), sep='\t', header=TRUE)

	# get SNPs into coloc format
	gfilt[,ref_snp := paste(chr, snp_pos, sep='_')]
	lfilt[,ref_snp := paste(chr, snp_pos, sep='_')]

	# all tests should overlap between global and local. first, merge global and local
	allp <- merge(gfilt, lfilt, by=c('gene','chr','snp_pos','ref_snp','ref','alt','build'), suffixes=c('.global','.local'))
	allp[,tissue := t]

	# now merge with coloc
	m <- merge(allp, coloc, by.x=c('gene','ref_snp', 'tissue'), by.y=c('feature','ref_snp','tissue'))

	m[,method := sapply(eqtl_file, get_method)]
	m[,eqtl_file := NULL]
	m[,base_gwas_file := NULL]

	# now I need to merge this back with the egenes file to see what the top SNPs were 
	egenes <- egenes_master[tissue==t]
	egenes <- egenes[pval_nominal < 1e-5]

	# filter egenes down to those in both methods with different lead SNPs 
	g <- egenes[method=='global']
	l <- egenes[method=='LAVA']

	g[,pair:=paste0(variant_id,':',gene_id)]
	l[,pair:=paste0(variant_id,':',gene_id)]
	repeats <- unique(l[pair %in% g[,pair], gene_id])

	g <- g[!(gene_id %in% repeats)]
	l <- l[!(gene_id %in% repeats)]

	egenes <- data.table(rbind(l, g))
	egenes[,ref_snp := sapply(variant_id, get_ref_snp)]
	egenes[,count:=NULL]
	egenes[,gtex_cutoff:=NULL]
	egenes[,slope:=NULL]

	master <- merge(egenes, m, by.x=c('gene_id','ref_snp','tissue','method'), by.y=c('gene','ref_snp','tissue','method'))

	master <- master[clpp>0.02]
	print(t)
	print(master)
	system('echo')
	#master[method.filt!=method.coloc] # are these top SNPs in global that are still significant in LAVA but just not lead SNPs? 

}
	


