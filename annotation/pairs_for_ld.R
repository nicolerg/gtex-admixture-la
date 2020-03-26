library(data.table)

# calculate LD from a pair of variants from ALL genes

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

# remove eGenes with overlapping lead SNP sets
tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

dflist <- list()
i <- 1 
for (t in tissues){

	egenes <- egenes_master[tissue==t]
	global <- egenes[method=='global']
	lava <- egenes[method=='LAVA']

	lava[,pair:=paste(gene_id, variant_id, sep=':')]
	global[,pair:=paste(gene_id, variant_id, sep=':')]
	repeats <- lava[pair %in% global[,pair], gene_id]

	egenes <- egenes[!gene_id %in% repeats]
	dflist[[i]] <- egenes

	i <- i + 1

}

egenes_filt <- rbindlist(dflist)

# concatenate variants
egenes <- egenes_filt[,list(lead_variants=paste(variant_id, collapse=';')),by=.(gene_id, count, tissue, method, pval_nominal)]

# get the first one from each tied set
g <- egenes[method == 'global']
l <- egenes[method == 'LAVA']
l[,method:=NULL]
g[,method:=NULL]

egenes = merge(g, l, by=c('gene_id','tissue','count'), suffixes=c('_global','_local'))
egenes[,lead_local := sapply(lead_variants_local, function(x) unname(unlist(strsplit(x,';')))[1])]
egenes[,lead_global := sapply(lead_variants_global, function(x) unname(unlist(strsplit(x,';')))[1])]
head(egenes)
nrow(egenes)

# write out a file 
write.table(egenes[,.(gene_id,tissue,lead_local,lead_global)], '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/20191030/evariant_pairs_test_ld.txt', sep='\t', col.names=F, row.names=F, quote=F)

#### RUN plink_ld.sh
