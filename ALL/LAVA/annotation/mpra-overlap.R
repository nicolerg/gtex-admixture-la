library(data.table)
library(ggplot2)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

# mpra <- '/users/nsabell/mpra/1KG/release/1KG_DESeq2_results.txt'
# mpra <- fread(mpra, sep='\t', header=TRUE)
# mpra <- mpra[,.(chr, pos, ref, alt, pvalue_expr)]
# head(mpra)

# map <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.rsid.map.txt.gz'
# map <- fread(cmd=sprintf('zcat %s',map), sep='\t', header=TRUE)
# head(map)

# map[, ref := sapply(hg38, function(x) unname(unlist(strsplit(x, '_')))[3])]
# map[, alt := sapply(hg38, function(x) unname(unlist(strsplit(x, '_')))[4])]

# mpra[, hg19 := paste0('chr',chr,'_',pos)]
# head(map)

# master <- merge(map, mpra, by=c('ref','alt','hg19'))
# head(master)
# nrow(master)

# save(master, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/mpra/merged.RData')
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/mpra/merged.RData')
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

# egenes_master <- egenes_master[pval_nominal < 1e-6]

dtlist <- list()
i<-1 
for (t in tissues){

	dt <- egenes_master[tissue==t]
	dt[,pair := paste(gene_id, variant_id, sep=':')]

	global <- dt[method=='global']
	lava <- dt[method=='LAVA']

	# remove eGenes that have overlapping top SNPs
	repeats <- unique(lava[pair %in% global[,pair],gene_id])
	global <- global[!(gene_id %in% repeats)]
	lava <- lava[!(gene_id %in% repeats)]

	# # remove eGenes that are not in both methods
	# g <- global[gene_id %in% lava[,gene_id]]
	# l <- lava[gene_id %in% global[,gene_id]]

	m <- data.table(rbind(global,lava))

	dtlist[[i]] <- m
	i <- i+1
}

egenes_master <- data.table(rbindlist(dtlist))

colnames(master)[colnames(master) == 'pvalue_expr'] <- 'mpra_pvalue_expr'
master <- master[,.(rsid, hg38, mpra_pvalue_expr)]
mpra_master <- merge(egenes_master, master, by.x='variant_id', by.y='hg38', all.x=TRUE)

subset <- mpra_master[,list(mpra_pvalue_expr=min(mpra_pvalue_expr, na.rm=TRUE)), by = .(gene_id,tissue,method,pval_nominal)]
table(subset[,method])
table(subset[mpra_pvalue_expr!=Inf, method])

t.test(subset[method=='global'&mpra_pvalue_expr!=Inf, mpra_pvalue_expr], subset[method=='LAVA'&mpra_pvalue_expr!=Inf, mpra_pvalue_expr])

# distributions of MPRA p-values

jpeg('~/gtex-admix/plots/final_figures/mpra-hist.jpg',res=150,width=1500,height=600)
ggplot(subset, aes(x=-log10(mpra_pvalue_expr), fill=factor(method))) +
	geom_histogram(alpha=0.3, bins=800, position='identity') +
	theme_bw() +
	labs(x='MPRA differential expression p-value (-log10)', y='Frequency') +
	scale_fill_manual(values=c(global='blue',LAVA='red'), name='Method', labels=c(global='Global',LAVA='Local')) +
	coord_cartesian(xlim=c(0,18))
dev.off()



