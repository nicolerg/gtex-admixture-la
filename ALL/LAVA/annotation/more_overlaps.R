library(data.table)
library(ggplot2)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

shortlab <- c(Adipose_Subcutaneous='Adipose',
	Artery_Tibial='Tibial\nartery',
	Lung='Lung',
	Muscle_Skeletal='Skeletal\nmuscle',
	Nerve_Tibial='Tibial\nnerve',
	Skin_Not_Sun_Exposed_Suprapubic='Skin (not\nsun-exposed)')

mcnemar <- function(a,b,c,d){
	m <- matrix(c(a,b,c,d), byrow=TRUE, nrow=2, ncol=2)
	p_val <- mcnemar.test(m)$p.value
	return(p_val)
}

# # Hunter Fraser caQTLs ====================================================================================================

caqtl <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/caQTL/elife-39595-supp2-v1.txt'
caqtl <- fread(caqtl, sep='\t', header=TRUE)

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

map <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/master/rsid_to_gtex_snp_map.tsv'
map <- fread(map, sep='\t', header=FALSE)
colnames(map) <- c('rsid', 'variant_id')

caqtl <- merge(caqtl, map, by='rsid')
master <- merge(egenes_master, caqtl, by='variant_id', all.x=TRUE)

caQTL <- function(cutoff){
	
	egenes <- master[pval_nominal < cutoff]

	# same eGenes, different lead SNPs 
	egenes[,pair := paste(gene_id, variant_id, sep=':')]

	dflist <- list()
	counts_df <- list()

	i <- 1
	for (t in tissues){

		dt <- egenes[tissue == t]

		global <- dt[method=='global']
		lava <- dt[method=='LAVA']

		# same eGenes, different lead SNPs 
		global <- global[gene_id %in% lava[,gene_id]]
		lava <- lava[gene_id %in% global[,gene_id]]
		repeats <- unique(lava[pair %in% global[,pair],gene_id])
		global <- global[!(gene_id %in% repeats)]
		lava <- lava[!(gene_id %in% repeats)]

		global[,annotation := ifelse(is.na(fishers_pval),0,1)]
		lava[,annotation := ifelse(is.na(fishers_pval),0,1)]

		global <- global[,.(gene_id,annotation)]
		lava <- lava[,.(gene_id,annotation)]
		
		g_annot <- global[,list(summ_annot=sum(annotation)), by = gene_id]
		g_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
		g_annot[,summ_annot:=NULL]
		l_annot <- lava[,list(summ_annot=sum(annotation)), by = gene_id]
		l_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
		l_annot[,summ_annot:=NULL]

		# merge 
		all <- merge(g_annot, l_annot, by='gene_id', all=TRUE, suffixes=c('.global','.lava'))

		counts_df[[i]] <- all

		total <- nrow(all)
		g_hit <- sum(all[,bin.global])
		l_hit <- sum(all[,bin.lava])

		df <- data.frame(tissue=t,
			method=c('global','LAVA'),
			hit=c(g_hit, l_hit),
			count=total)

		dflist[[i]] <- df

		i <- i + 1

	}

	m <- data.table(rbindlist(dflist))
	labels <- m[,.(tissue, count)]
	labels <- unique(labels)
	print(labels)

	label_pos <- 1.1*max(m[,hit])

	g <- ggplot() +
		geom_bar(data=m, aes(x=factor(tissue),y=hit,fill=factor(method)), stat='identity', position='dodge') +
		geom_text(data=labels, aes(label = paste0('N = ',count), x = tissue, y = label_pos), vjust = -0.5) +
		scale_fill_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='Local'), name='Method') +
		labs(y='N top SNP sets overlapping caQTLs',title=cutoff) +
		theme_bw() + 
		theme(axis.title.x=element_blank(),
			axis.text.x=element_text(colour='black',size=10)) +
		scale_x_discrete(labels=shortlab) +
		ylim(0,1.2*max(m[,hit])) 

	pdf(paste0("~/gtex-admix/plots/final_figures/caQTL-same-",cutoff,".pdf"),width=10,height=4)
	print(g)
	dev.off()

	# get pval across all tissues 
	counts <- data.table(rbindlist(counts_df))
	counts <- unique(counts)

	tab <- table(counts[,bin.global], counts[,bin.lava])
	print(tab)
	a <- tab[1,1]
	b <- tab[1,2]
	c <- tab[2,1]
	d <- tab[2,2]
	pval <- mcnemar(a,b,c,d)
	print(pval)

}

# caQTL(1e-6)
# caQTL(1e-5)

# Hunter Fraser bQTLs ====================================================================================================

bqtl <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/bQTL/bQTL.merged.txt'
bqtl <- fread(bqtl, sep='\t', header=TRUE)

# convert to hg38
build_map <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/master/gtex_hg19/egenes_master.hg19.to.hg39.map.txt',sep='\t',header=TRUE)

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')
egenes <- merge(egenes_master, build_map, by.x='variant_id', by.y='hg38')

bqtl[,hg19 := paste(Chr, position, sep='_')]

master <- merge(egenes, bqtl, by='hg19', all.x=TRUE)
head(master[complete.cases(master)])

master[,annotation := ifelse(is.na(pvalue),0,1)]

bQTL <- function(cutoff){
	
	egenes <- master[pval_nominal < cutoff]

	# same eGenes, different lead SNPs 
	egenes[,pair := paste(gene_id, variant_id, sep=':')]

	dflist <- list()
	counts_df <- list()

	i <- 1
	for (t in tissues){

		dt <- egenes[tissue == t]

		global <- dt[method=='global']
		lava <- dt[method=='LAVA']

		# same eGenes, different lead SNPs 
		global <- global[gene_id %in% lava[,gene_id]]
		lava <- lava[gene_id %in% global[,gene_id]]
		repeats <- unique(lava[pair %in% global[,pair],gene_id])
		global <- global[!(gene_id %in% repeats)]
		lava <- lava[!(gene_id %in% repeats)]

		global <- global[,.(gene_id,annotation)]
		lava <- lava[,.(gene_id,annotation)]
		
		g_annot <- global[,list(summ_annot=sum(annotation)), by = gene_id]
		g_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
		g_annot[,summ_annot:=NULL]
		l_annot <- lava[,list(summ_annot=sum(annotation)), by = gene_id]
		l_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
		l_annot[,summ_annot:=NULL]

		# merge 
		all <- merge(g_annot, l_annot, by='gene_id', all=TRUE, suffixes=c('.global','.lava'))

		counts_df[[i]] <- all

		total <- nrow(all)
		g_hit <- sum(all[,bin.global])
		l_hit <- sum(all[,bin.lava])

		df <- data.frame(tissue=t,
			method=c('global','LAVA'),
			hit=c(g_hit, l_hit),
			count=total)

		dflist[[i]] <- df

		i <- i + 1

	}

	m <- data.table(rbindlist(dflist))
	labels <- m[,.(tissue, count)]
	labels <- unique(labels)
	print(labels)

	label_pos <- 1.1*max(m[,hit])

	g <- ggplot() +
		geom_bar(data=m, aes(x=factor(tissue),y=hit,fill=factor(method)), stat='identity', position='dodge') +
		geom_text(data=labels, aes(label = paste0('N = ',count), x = tissue, y = label_pos), vjust = -0.5) +
		scale_fill_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='Local'), name='Method') +
		labs(y='N top SNP sets overlapping bQTLs',title=cutoff) +
		theme_bw() + 
		theme(axis.title.x=element_blank(),
			axis.text.x=element_text(colour='black',size=10)) +
		scale_x_discrete(labels=shortlab) +
		ylim(0,1.2*max(m[,hit])) 

	pdf(paste0("~/gtex-admix/plots/final_figures/bQTL-same-",cutoff,".pdf"),width=10,height=4)
	print(g)
	dev.off()

	# get pval across all tissues 
	counts <- data.table(rbindlist(counts_df))
	counts <- unique(counts)

	tab <- table(counts[,bin.global], counts[,bin.lava])
	if (nrow(tab) > 1){
		print(tab)
		a <- tab[1,1]
		b <- tab[1,2]
		c <- tab[2,1]
		d <- tab[2,2]
		pval <- mcnemar(a,b,c,d)
		print(pval)
	}


}

bQTL(1e-6)
bQTL(1e-5)
