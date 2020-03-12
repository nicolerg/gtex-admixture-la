# GTEx v8 local ancestry final paper figures

library(data.table)
library(ggplot2)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

tissuecols <- c(Adipose_Subcutaneous="#FF6600",
	Artery_Tibial="#FF0000",
	Lung="#99FF00",
	Muscle_Skeletal="#AAAAFF",
	Nerve_Tibial="#FFD700",
	Skin_Not_Sun_Exposed_Suprapubic="#0000FF")

tissuelabs <- c(Adipose_Subcutaneous="Adipose (subcutaneous)",
	Artery_Tibial="Tibial artery",
	Lung="Lung",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Skin_Not_Sun_Exposed_Suprapubic="Skin (not sun-exposed)")

shortlab <- c(Adipose_Subcutaneous='Adipose',
	Artery_Tibial='Tibial\nartery',
	Lung='Lung',
	Muscle_Skeletal='Skeletal\nmuscle',
	Nerve_Tibial='Tibial\nnerve',
	Skin_Not_Sun_Exposed_Suprapubic='Skin (not\nsun-exposed)')

chi_sq <- function(lava_hit, global_hit, total_lava, total_global){
	a <- total_lava - lava_hit
	b <- total_global - global_hit
	c <- lava_hit
	d <- global_hit 

	m <- matrix(c(a,b,c,d), byrow=TRUE, nrow=2, ncol=2)
	p_val <- chisq.test(m)$p.value
	return(p_val)
}

fisher <- function(lava_hit, global_hit, total_lava, total_global){
	a <- total_lava - lava_hit
	b <- total_global - global_hit
	c <- lava_hit
	d <- global_hit 

	m <- matrix(c(a,b,c,d), byrow=TRUE, nrow=2, ncol=2)
	p_val <- fisher.test(m)$p.value
	return(p_val)
}

mcnemar <- function(a,b,c,d){
	m <- matrix(c(a,b,c,d), byrow=TRUE, nrow=2, ncol=2)
	p_val <- mcnemar.test(m)$p.value
	return(p_val)
}

t_test <- function(master,val){

	ttest <- data.table(tissue=tissues,
		pvalue=1.0,
		greater_lava=1.0,
		greater_global=1.0)

	for (t in tissues){
		g <- master[tissue == t & method == 'global',get(val)]
		l <- master[tissue == t & method == 'LAVA',get(val)]
		g_global <- t.test(g,l,na.action=na.omit(),alternative='greater')$p.value
		g_lava <- t.test(g,l,na.action=na.omit(),alternative='less')$p.value

		pval <- t.test(g,l,na.action=na.omit())$p.value
		ttest[tissue==t, pvalue := pval]
		ttest[tissue==t, greater_lava := g_lava]
		ttest[tissue==t, greater_global := g_global]
	}

	return(ttest)
}

# allp_dt <- list()
# egenes_dt <- list()

# i <- 1

# for (t in tissues){

# 	basedir = paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/',t,'/')
# 	filt_allpairs = paste0(basedir,t,'-LAVA-global-merged-allpairs-filt.tsv.gz')
# 	global_egenes = paste0(basedir,'gtex.admix.global.egenes.tied.txt.gz')
# 	lava_egenes = paste0(basedir,'gtex.admix.lava.egenes.tied.txt.gz')

# 	allp_filt = fread(cmd=sprintf("zcat %s",filt_allpairs), sep='\t', header=TRUE)
# 	global = fread(cmd=sprintf("zcat %s",global_egenes), sep='\t', header=TRUE)
# 	lava = fread(cmd=sprintf("zcat %s",lava_egenes), sep='\t', header=TRUE)

# 	# make a merged allpairs file
# 	allp_filt[,tissue:=t]

# 	# make a merged egenes file 
# 	global[,method:='global']
# 	global[,tissue:=t]
# 	lava[,method:='LAVA']
# 	lava[,tissue:=t]
# 	egenes = data.table(rbind(global, lava))

# 	allp_dt[[i]] <- allp_filt
# 	egenes_dt[[i]] <- egenes

# 	i <- i + 1

# }

# allp_master <- data.table(rbindlist(allp_dt))
# egenes_master <- data.table(rbindlist(egenes_dt))
# egenes_master[,pval_nominal := as.numeric(pval_nominal)]

# save(allp_master, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/allp_master.RData')
# save(egenes_master, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/allp_master.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

# # eGene discovery, global vs LAVA (one point per tissue) ====================================================================
# Try different cutoffs: 1e-6, 1e-5

# df <- data.table(tissue=c(tissues,tissues), global=0, LAVA=0, pval_cutoff=c(rep(1e-5, length(tissues)), rep(1e-6, length(tissues))))

# for (t in tissues){
# 	for (cutoff in c(1e-5,1e-6)){

# 		sub <- egenes_master[ tissue == t & pval_nominal < cutoff ]
# 		df[tissue == t & pval_cutoff == cutoff, global := length(unique(sub[method == 'global', gene_id]))]
# 		df[tissue == t & pval_cutoff == cutoff, LAVA := length(unique(sub[method == 'LAVA', gene_id]))]

# 	}
# }

# save(df, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egene_discovery_df.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egene_discovery_df.RData')

# pdf("~/gtex-admix/plots/final_figures/egene_discovery.pdf",width=8,height=3)
# ggplot(df, aes(x=global, y=LAVA, fill=factor(tissue))) +
# 	geom_point(colour='black', shape=21) +
# 	geom_abline(slope=1, intercept=0, linetype='dashed') +
# 	scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue') +
# 	labs(x='N global eGenes', y='N LAVA eGenes') +
# 	facet_wrap(~ pval_cutoff, ncol=2, scales='free') +
# 	theme_bw() +
#	scale_x_discrete(labels=shortlab)
# dev.off()

# sub <- df[pval_cutoff == 1e-6]
# egene_discovery_1e6 <- ggplot(sub, aes(x=global, y=LAVA, fill=factor(tissue))) +
# 	geom_point(colour='black', shape=21) +
# 	geom_abline(slope=1, intercept=0, linetype='dashed') +
# 	scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue') +
# 	labs(x='N global eGenes', y='N LAVA eGenes') +
# 	theme_bw() +
#	scale_x_discrete(labels=shortlab)

# pdf("~/gtex-admix/plots/final_figures/egene_discovery_1e-6.pdf",width=6,height=3)
# egene_discovery_1e6
# dev.off()

# ## qq plot ========================================================================================================================================

# # one line per tissue 

# dtlist <- list()
# i <- 1
# for (t in tissues){
# 	sub <- allp_master[tissue == t]
# 	d <- as.data.table(qqplot(sub[,-log10(global_pval_nominal)], sub[,-log10(lava_pval_nominal)], plot.it=FALSE))
# 	colnames(d) <- c('global','lava')
# 	d[,tissue:=t]
# 	dtlist[[i]] <- d 
# 	i <- i+1 
# }

# # qq <- data.table(rbindlist(dtlist))

# # save(qq, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/qq.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/qq.RData')

# jpeg("~/gtex-admix/plots/final_figures/qq_plot.jpg",width=1500,height=900,res=200)
# ggplot(qq, aes(x=global, y=lava, colour=factor(tissue))) +
# 	geom_line() +
# 	scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue') +
# 	labs(x='Global nominal p-value (-log10)', y='LAVA nominal p-value (-log10)') +
# 	geom_abline(slope=1, intercept=0, linetype='dashed')  +
# 	theme_bw() +
#	scale_x_discrete(labels=shortlab)
# dev.off()

# # TSS distance ========================================================================================================================================

# plot TSS distances 

# first, 1e-5

# x <- allp_master[, .(gene_id,variant_id,tss_distance) ]
# tss <- merge(egenes_master, x, by=c('gene_id','variant_id'))
# head(tss)

# save(tss,file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss.RData')

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss.RData')

# # first, for 1e-5, only matching eGenes

# tss <- tss[!duplicated(tss)] # not sure why so many rows are duplicated 
# tss[,pval_nominal := as.numeric(pval_nominal)]

# tss_same_egenes <- function(pval_cutoff){
# 	tss_sub <- tss[pval_nominal < pval_cutoff]

# 	tss_df <- list()
# 	i <- 1

# 	for (t in tissues){

# 		tss_min <- data.table(gene_id=unique(tss_sub[,gene_id]),lava=1e7,global=1e7)
# 		sub1 <- tss_sub[tissue == t]
# 		print(table(sub1[,method]))

# 		for (g in unique(sub1[,gene_id])){
			
# 			lava <- sub1[gene_id==g & method=='LAVA']
# 			global <- sub1[gene_id==g & method=='global']

# 			# skip if gene not present in both LAVA and global
# 			if(nrow(lava) == 0){ next }
# 			lava_min <- lava[ which.min(abs(lava[,tss_distance])), tss_distance ]
# 			if(nrow(global) == 0){ next }
# 			global_min <- global[ which.min(abs(global[,tss_distance])), tss_distance ]

# 			# skip if same variant 
# 			if(global_min == lava_min){ next } 

# 			tss_min[gene_id == g, lava := lava_min]
# 			tss_min[gene_id == g, global := global_min]
# 		}

# 		tss_min[,tissue:=t]
# 		tss_df[[i]] <- tss_min
# 		i <- i + 1 
# 	}

# 	tss_master <- data.table(rbindlist(tss_df))
# 	tss_master[tss_master == 1e7] <- NA

# 	return(tss_master)
# }

# tss_uniq_too <- function(pval_cutoff){
# # try a version where unique egenes are included 

# 	tss_sub <- tss[pval_nominal < pval_cutoff]

# 	tss_df <- list()
# 	i <- 1

# 	for (t in tissues){

# 		tss_min <- data.table(gene_id=unique(tss_sub[,gene_id]),lava=1e7,global=1e7)
# 		sub1 <- tss_sub[tissue == t]
# 		print(table(sub1[,method]))

# 		for (g in unique(sub1[,gene_id])){
			
# 			lava <- sub1[gene_id==g & method=='LAVA']
# 			global <- sub1[gene_id==g & method=='global']

# 			# keep unique genes
# 			if(nrow(lava) == 0){ 
# 				lava_min <- 1e7
# 			} else {
# 				lava_min <- lava[ which.min(abs(lava[,tss_distance])), tss_distance ]
# 			}
			
# 			if(nrow(global) == 0){ 
# 				global_min <- 1e7
# 			} else {
# 				global_min <- global[ which.min(abs(global[,tss_distance])), tss_distance ]
# 			}
			
# 			# skip if same variant 
# 			if(global_min == lava_min){
# 				next 
# 			}
			 
# 			tss_min[gene_id == g, lava := lava_min]
# 			tss_min[gene_id == g, global := global_min]
# 		}

# 		tss_min[,tissue:=t]
# 		tss_df[[i]] <- tss_min
# 		i <- i + 1 
# 	}

# 	tss_master_plus_uniq <- data.table(rbindlist(tss_df))
# 	tss_master_plus_uniq[tss_master_plus_uniq == 1e7] <- NA 
# 	return(tss_master_plus_uniq)
# }

# tss_master_1e5 <- tss_same_egenes(1e-5)
# save(tss_master_1e5, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_1e5.RData')

# tss_master_plus_uniq_1e5 <- tss_uniq_too(1e-5)
# save(tss_master_plus_uniq_1e5, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_plus_uniq_1e5.RData')

# # now for 1e-6

# tss_master_1e6 <- tss_same_egenes(1e-6)
# save(tss_master_1e6, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_1e6.RData')

# tss_master_plus_uniq_1e6 <- tss_uniq_too(1e-6)
# save(tss_master_plus_uniq_1e6, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_plus_uniq_1e6.RData')

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_1e5.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_plus_uniq_1e5.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_1e6.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/tss_master_plus_uniq_1e6.RData')

# # now plot

# for (df in c('tss_master_plus_uniq_1e5','tss_master_plus_uniq_1e6')){
# 	print(df)
# 	dt <- get(df)

# 	dt <- melt(dt, id=c('gene_id','tissue'))
# 	dt <- dt[complete.cases(dt)]
# 	counts <- as.data.table(table(dt[,variable]))
# 	print(as.data.table(table(dt[,variable],dt[,tissue])))
# 	lava_count <- counts[V1=='lava',N][1]
# 	global_count <- counts[V1=='global',N][1]

# 	g1 <- ggplot(dt, aes(x=value, fill=factor(variable))) +
# 		geom_density(alpha=0.2) +
# 		scale_fill_manual(values=c(global='blue',lava='red'),name='Method',
# 			labels=c(global=paste0('Global (N=',global_count,')'),lava=paste0('LAVA (N=',lava_count,')'))) +
# 		labs(title='All eGenes with different top SNP',x='TSS distance',y='Density') +
# 		geom_vline(xintercept=0,linetype="dashed") +
# 		facet_wrap(~tissue,ncol=3)
	
# 	pdf(paste0("~/gtex-admix/plots/final_figures/",df,".pdf"),width=10,height=6)
# 	print(g1)
# 	dev.off()
# }

# for (df in c('tss_master_1e5','tss_master_1e6')){
# 	print(df)
# 	dt <- get(df)

# 	dt <- melt(dt, id=c('gene_id','tissue'))
# 	dt <- dt[complete.cases(dt)]
# 	counts <- as.data.table(table(dt[,variable]))
# 	print(as.data.table(table(dt[,variable],dt[,tissue])))
# 	lava_count <- counts[V1=='lava',N][1]
# 	global_count <- counts[V1=='global',N][1]

# 	g1 <- ggplot(dt, aes(x=value, fill=factor(variable))) +
# 		geom_density(alpha=0.2) +
# 		scale_fill_manual(values=c(global='blue',lava='red'),name='Method',
# 			labels=c(global=paste0('Global (N=',global_count,')'),lava=paste0('LAVA (N=',lava_count,')'))) +
# 		labs(title='Same eGenes, different top SNPs',x='TSS distance',y='Density') +
# 		geom_vline(xintercept=0,linetype="dashed") +
# 		facet_wrap(~tissue,ncol=3)
	
# 	pdf(paste0("~/gtex-admix/plots/final_figures/",df,".pdf"),width=10,height=6)
# 	print(g1)
# 	dev.off()
# }

# ## GWAS coloc ========================================================================================================================================

# # generated by gwas_coloc.py

coloc_dt <- list()
i <- 1
for (t in tissues){
	file <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_egene_to_gwas_map_',t,'.tsv')
	dt <- fread(file, sep='\t', header=TRUE)
	coloc_dt[[i]] <- dt
	i <- i+1
}
coloc_master <- data.table(rbindlist(coloc_dt))
coloc_master <- melt(coloc_master, id=c('gene_id','tissue'))

coloc_plot <- function(cutoff){

	# egenes_per_tissue <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/master/uniq_egenes_per_tissue_',cutoff,'.tsv')

	# egenes_uniq <- fread(egenes_per_tissue, sep='\t', header=TRUE)

	# coloc <- merge(egenes_uniq, coloc_master, by=c('gene_id','tissue'), all.x=TRUE)
	# coloc[,value := -log10(value)]

	# # plot all diseases per tissue

	# # # get t.test p-values for each tissue
	# # coloc_sig <- t_test(coloc,'value')

	# # g1 <- ggplot(coloc_sig, aes(x=factor(tissue), y=value)) +
	# # 	geom_boxplot(aes(colour=factor(method))) +
	# # 	geom_text(aes(label = signif(pvalue,digits=3), x = tissue, y = 175)) +
	# # 	labs(y='GWAS colocalization p-value (-log10)') + 
	# # 	scale_colour_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
	# # 	theme(axis.title.x = element_blank()) +
	# # 	scale_x_discrete(labels=shortlab)

	# # pdf(paste0("~/gtex-admix/plots/final_figures/gwas-coloc-all-pheno-",cutoff,"-uniq-egenes.pdf"),width=10,height=6)
	# # print(g1)
	# # dev.off()

	# # # split out by disease

	# # g2 <- ggplot(coloc, aes(x=factor(variable), y=value, colour=factor(method))) +
	# # 	geom_boxplot() +
	# # 	labs(y='GWAS colocalization p-value (-log10)') + 
	# # 	scale_colour_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
	# # 	theme(axis.title.x = element_blank(),
	# # 		axis.text.x = element_text(angle=90,hjust=1)) +
	# # 	facet_wrap(~tissue, ncol=3,nrow=2)

	# # pdf(paste0("~/gtex-admix/plots/final_figures/gwas-coloc-bytissue-",cutoff,"-1e-5-uniq-egenes.pdf"),width=30,height=20)
	# # print(g2)
	# # dev.off()


	# # get t.test p-values for each disease, each tissue

	# ttest <- data.table(tissue=rep(tissues, length(unique(coloc[,variable]))),
	# 	phenotype=rep(unique(coloc[,variable]), each=length(tissues)),
	# 	pvalue=1.0,
	# 	greater_lava=1.0,
	# 	greater_global=1.0)

	# for (t in tissues){
	# 	tsub <- coloc[tissue == t]
	# 	phenos <- unique(tsub[,variable])
	# 	for (p in phenos){
	# 		g <- tsub[variable == p & method == 'global',value]
	# 		l <- tsub[variable == p & method == 'LAVA',value]
	# 		g_global <- t.test(g,l,na.action=na.omit(),alternative='greater')$p.value
	# 		g_lava <- t.test(g,l,na.action=na.omit(),alternative='less')$p.value

	# 		pval <- t.test(g,l,na.action=na.omit())$p.value
	# 		ttest[tissue==t&phenotype==p, pvalue := pval]
	# 		ttest[tissue==t&phenotype==p, greater_lava := g_lava]
	# 		ttest[tissue==t&phenotype==p, greater_global := g_global]
	# 	}
	# }

	# save(ttest, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/coloc_ttest.RData')
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/coloc_ttest.RData')

	print(head(ttest))

	print(ttest[pvalue < 0.05])

	ttest[,pvalue := -log10(pvalue)]
	ttest[,pvalue := ifelse(greater_global < greater_lava, -1*pvalue, pvalue)]
	ttest <- ttest[order(pvalue, decreasing=TRUE)]
	ttest[,index := 1:nrow(ttest)]
	ttest[,colour := ifelse(pvalue > 0, 'red', 'blue')]

	pdf('~/gtex-admix/plots/final_figures/gwas-coloc-p4.pdf', width=8, height=3)
	g<-ggplot(ttest, aes(x=index, y=pvalue, fill=colour)) +
		geom_bar(stat='identity') + 
		theme_bw() +
		theme(axis.text.x = element_blank(),
			axis.title.x = element_blank(),
			axis.title.y = element_blank()) +
		geom_hline(yintercept=-log10(0.05), linetype='dashed') +
		geom_hline(yintercept=log10(0.05), linetype='dashed') +
		geom_hline(yintercept=0) +
		scale_fill_manual(values=c(red='red',blue='blue'), guide='none')
	print(g)
	dev.off()

	nrow(ttest)
	npheno <- length(unique(coloc[,variable]))
	mt_pval <- 0.05/npheno

	print(mt_pval)
	ttest <- ttest[pvalue < 0.05]

	filt_coloc <- merge(coloc, ttest, by.x=c('tissue','variable'), by.y=c('tissue','phenotype'))

	filt_coloc[,label := paste(tissue,variable,sep=' : ')]

	g3 <- ggplot(filt_coloc, aes(x=factor(label), y=value)) +
		geom_boxplot(aes(colour=factor(method))) +
		geom_text(aes(label = signif(pvalue,digits=3), x = label, y = 80)) +
		labs(y='GWAS colocalization p-value (-log10)') + 
		scale_colour_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
		theme(axis.title.y = element_blank()) +
		coord_flip()

	pdf(paste0("~/gtex-admix/plots/final_figures/gwas-filt-coloc-",cutoff,"-uniq-egenes.pdf"),width=12,height=8)
	print(g3)
	dev.off()

	print(length(unique(coloc[,variable])))
}

# coloc_plot(1e-5)
# coloc_plot(1e-6)


# ## CADD ========================================================================================================================================

# generated by parse_cadd.py

cadd <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cadd/gtex_snp_to_cadd_map.tsv'
cadd <- fread(cadd, sep='\t', header=FALSE)
colnames(cadd) <- c('variant_id','cadd_score')
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

cadd_plot <- function(cutoff){

	print(cutoff)

	# read in unique eQTLs (includes diff lead SNP)
	load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/uniq_eqtl_',cutoff,'.RData'))
	cadd_master <- merge(master_uniq_eqtl, cadd, by='variant_id', all.x=TRUE)

	c <- cadd_master[,.(gene_id, cadd_score, tissue, method)]
	c <- c[,list(max_score=max(cadd_score)), by = .(gene_id,tissue,method)]

	# get t.test p-values for each tissue
	cadd_sig <- t_test(c,'max_score')

	labels <- cadd_sig[,.(tissue, pvalue)]
	labels <- unique(labels)
	labels[,y := max(c[,max_score],na.rm=TRUE)+0.1]

	print('All unique eQTLs')
	print(t.test(c[method=='global',max_score],c[method=='LAVA',max_score])$p.value)

	g1 <- ggplot() +
		geom_boxplot(data=c, aes(x=factor(tissue), y=max_score, colour=factor(method))) +
		geom_text(data=labels, aes(x=factor(tissue), y=y, label=signif(pvalue,3)),col='black') +
		labs(y='CADD score') + 
		scale_colour_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
		theme_bw() + 
		theme(axis.title.x = element_blank(),
			axis.text.x = element_text(colour='black')) +
		scale_x_discrete(labels=shortlab)


	pdf(paste0("~/gtex-admix/plots/final_figures/cadd-score-",cutoff,"-all-egenes.pdf"),width=6,height=4)
	print(g1)
	dev.off()

	# now for same genes where lead SNP differs (master_diff_lead)

	load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_',cutoff,'.RData'))
	cadd_master <- merge(master_diff_lead, cadd, by='variant_id', all.x=TRUE)

	c <- cadd_master[,.(gene_id, cadd_score, tissue, method)]
	c <- c[,list(max_score=max(cadd_score)), by = .(gene_id,tissue,method)]

	# get t.test p-values for each tissue
	cadd_sig <- t_test(c,'max_score')

	labels <- cadd_sig[,.(tissue, pvalue)]
	labels <- unique(labels)
	labels[,y := max(c[,max_score],na.rm=TRUE)+0.1]

	print('Same eGene, different lead SNP')
	print(t.test(c[method=='global',max_score],c[method=='LAVA',max_score])$p.value)

	g2 <- ggplot() +
		geom_boxplot(data=c, aes(x=factor(tissue), y=max_score, colour=factor(method))) +
		geom_text(data=labels, aes(x=factor(tissue), y=y, label=signif(pvalue,3)),col='black') +
		labs(y='CADD score') + 
		theme_bw() +
		scale_colour_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
		theme(axis.title.x = element_blank(),
			axis.text.x = element_text(colour='black')) +
		scale_x_discrete(labels=shortlab)

	pdf(paste0("~/gtex-admix/plots/final_figures/cadd-score-",cutoff,"-same-egenes.pdf"),width=6,height=4)
	print(g2)
	dev.off()
}

cadd_plot(1e-5)
cadd_plot(1e-6)
q()


# # Fst ========================================================================================================================================

# # output of gtex_v8_rsid_map.py, fst.py 
# fst <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/gtex_snp_to_fst_map.tsv'
# fst_map <- fread(fst, header=TRUE, sep='\t')
# # columns in published data are mislabelled 
# colnames(fst_map) <- c('variant_id','rsid','fst_score', 'fst_pvalue')

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

# fst <- function(cutoff){

# 	egenes_m <- egenes_master[pval_nominal < cutoff]
# 	merge_fst <- merge(egenes_m, fst_map, by='variant_id', all.x=TRUE)

# 	# # first, for all egenes with different lead SNPs

# 	# dtlist <- list()
# 	# countdf <- list()
# 	# i<-1 
# 	# for (t in tissues){

# 	# 	dt <- merge_fst[tissue==t]
# 	# 	dt[,pair := paste(gene_id, variant_id, sep=':')]

# 	# 	global <- dt[method=='global']
# 	# 	lava <- dt[method=='LAVA']

# 	# 	# remove eGenes that have overlapping top SNPs
# 	# 	repeats <- unique(lava[pair %in% global[,pair],gene_id])
# 	# 	global <- global[!(gene_id %in% repeats)]
# 	# 	lava <- lava[!(gene_id %in% repeats)]

# 	# 	# need to keep track of how many values there are for each tissue and method 
# 	# 	m <- data.table(rbind(global,lava))

# 	# 	m[,annotation := ifelse(is.na(fst_score),0,1)]

# 	# 	dtlist[[i]] <- m
# 	# 	i <- i+1
# 	# }

# 	# master_fst <- data.table(rbindlist(dtlist))
# 	# master_fst <- master_fst[,.(gene_id, fst_score, tissue, method)]
# 	# master_fst <- master_fst[,list(max_score=max(fst_score, na.rm=TRUE)), by = .(gene_id,tissue,method)]
# 	# master_fst[master_fst == -Inf] <- NA
# 	# master_fst <- master_fst[complete.cases(master_fst)]
	
# 	# counts <- as.data.frame.matrix(table(master_fst[,tissue],master_fst[,method]))
# 	# counts <- data.table(counts,keep.rownames=TRUE)
# 	# counts <- melt(counts, id='rn')

# 	# pdf(paste0("~/gtex-admix/plots/final_figures/fst-",cutoff,"-all.pdf"),width=10,height=6)
# 	# g<-ggplot() +
# 	# 	geom_boxplot(data=master_fst,aes(x=factor(tissue),y=max_score,colour=factor(method))) +
# 	# 	geom_text(data=counts, aes(x=factor(rn), y=1, label=value, colour=factor(variable)), position = position_dodge(width = 0.9)) +
# 	# 	scale_colour_manual(values=c(LAVA='red', global='blue'),
# 	# 		labels=c(LAVA='LAVA',global='Global'), name='Method') +
# 	# 	labs(y='Fst score') +
# 	# 	theme(axis.title.x=element_blank()) +
# 	#	scale_x_discrete(labels=shortlab)
# 	# print(g)
# 	# dev.off()

# 	# # take fst_pvalue into consideration?

# 	# master_fst <- master_fst[fst_pvalue < 0.05]
# 	# counts <- as.data.frame.matrix(table(master_fst[annotation==1,tissue],master_fst[annotation==1,method]))
# 	# counts <- data.table(counts,keep.rownames=TRUE)
# 	# counts <- melt(counts, id='rn')

# 	# pdf(paste0("~/gtex-admix/plots/final_figures/fst-",cutoff,"-p05-all.pdf"),width=10,height=6)
# 	# g<-ggplot() +
# 	# 	geom_boxplot(data=master_fst,aes(x=factor(tissue),y=fst_score,colour=factor(method))) +
# 	# 	geom_text(data=counts, aes(x=factor(rn), y=0.6, label=value, colour=factor(variable)), position = position_dodge(width = 0.9)) +
# 	# 	scale_colour_manual(values=c(LAVA='red', global='blue'),
# 	# 		labels=c(LAVA='LAVA',global='Global'), name='Method') +
# 	# 	labs(y='Fst score') +
# 	# 	theme(axis.title.x=element_blank()) +
# 	#	scale_x_discrete(labels=shortlab)
# 	# 	ylim(0,0.6)
# 	# print(g)
# 	# dev.off()


# 	# next, for matched egenes with different lead SNPs ==============================

# 	dtlist <- list()
# 	countdf <- list()
# 	i<-1 
# 	for (t in tissues){

# 		dt <- merge_fst[tissue==t]
# 		dt[,pair := paste(gene_id, variant_id, sep=':')]
# 		dt <- dt[complete.cases(dt)]

# 		global <- dt[method=='global']
# 		lava <- dt[method=='LAVA']

# 		# remove eGenes that have overlapping top SNPs
# 		repeats <- unique(lava[pair %in% global[,pair],gene_id])
# 		global <- global[!(gene_id %in% repeats)]
# 		lava <- lava[!(gene_id %in% repeats)]

# 		# remove eGenes that are not in both methods
# 		global <- global[gene_id %in% lava[,gene_id]]
# 		lava <- lava[gene_id %in% global[,gene_id]]

# 		# keep only those whose difference in p-value is above some cutoff 

# 		global <- global[,.(gene_id, pval_nominal, tissue, fst_score)]
# 		lava <- lava[,.(gene_id, pval_nominal, tissue, fst_score)]

# 		global <- global[,list(max_score=fst_score[which.max(abs(fst_score))]), by = .(gene_id,tissue,pval_nominal)]
# 		lava <- lava[,list(max_score=fst_score[which.max(abs(fst_score))]), by = .(gene_id,tissue,pval_nominal)]

# 		m <- merge(global, lava, by=c('gene_id','tissue'), suffixes=c('.global','.lava'))
# 		m_filt <- m[abs(-log10(pval_nominal.global) - -log10(pval_nominal.lava))>=1]

# 		dtlist[[i]] <- m_filt
# 		i <- i+1

# 		print(t.test(m_filt[,max_score.lava], m_filt[,max_score.global]))

# 	}

# 	master_fst <- data.table(rbindlist(dtlist))
# 	counts <- as.data.table(table(master_fst[,tissue]))
# 	colnames(counts) <- c('tissue', 'N')

# 	print(t.test(master_fst[,max_score.lava], master_fst[,max_score.global]))

# 	master_fst <- melt(master_fst, measure.vars=c('max_score.global', 'max_score.lava'))

# 	save(master_fst, counts, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/figure_data/master_fst.RData')

# 	pdf(paste0("~/gtex-admix/plots/final_figures/fst-",cutoff,"-same.pdf"),width=8,height=4)
# 	g<-ggplot() +
# 		geom_boxplot(data=master_fst,aes(x=factor(tissue),y=value,fill=factor(variable))) +
# 		geom_text(data=counts, aes(x=factor(tissue), y=1, label=paste0('N = ',N)), colour='black') +
# 		scale_fill_manual(values = c('max_score.global'='blue', 'max_score.lava'='red'),
# 			labels=c('max_score.lava'='Local','max_score.global'='Global'), name='Method') +
# 		labs(y='Fst score') +
# 		theme_bw() +
# 		theme(axis.title.x=element_blank(),
# 			axis.text.x=element_text(colour='black')) +
# 		scale_x_discrete(labels=shortlab)
# 	print(g)
# 	dev.off()

# 	# pdf(paste0("~/gtex-admix/plots/final_figures/fst-",cutoff,"-same.pdf"),width=10,height=6)
# 	# g<-ggplot() +
# 	# 	geom_boxplot(data=master_fst,aes(x=factor(tissue),y=fst_score,colour=factor(method))) +
# 	# 	geom_text(data=counts, aes(x=factor(rn), y=1, label=value, colour=factor(variable)), position = position_dodge(width = 0.9)) +
# 	# 	scale_colour_manual(values=c(LAVA='red', global='blue'),
# 	# 		labels=c(LAVA='LAVA',global='Global'), name='Method') +
# 	# 	labs(y='Fst score') +
# 	# 	theme(axis.title.x=element_blank()) +
# 	#	scale_x_discrete(labels=shortlab)
# 	# print(g)
# 	# dev.off()

# 	# # take fst_pvalue into consideration?

# 	# master_fst <- master_fst[fst_pvalue < 0.05]
# 	# counts <- as.data.frame.matrix(table(master_fst[annotation==1,tissue],master_fst[annotation==1,method]))
# 	# counts <- data.table(counts,keep.rownames=TRUE)
# 	# counts <- melt(counts, id='rn')

# 	# pdf(paste0("~/gtex-admix/plots/final_figures/fst-",cutoff,"-p05-same.pdf"),width=10,height=6)
# 	# g<-ggplot() +
# 	# 	geom_boxplot(data=master_fst,aes(x=factor(tissue),y=fst_score,colour=factor(method))) +
# 	# 	geom_text(data=counts, aes(x=factor(rn), y=0.6, label=value, colour=factor(variable)), position = position_dodge(width = 0.9)) +
# 	# 	scale_colour_manual(values=c(LAVA='red', global='blue'),
# 	# 		labels=c(LAVA='LAVA',global='Global'), name='Method') +
# 	# 	labs(y='Fst score') +
# 	# 	theme(axis.title.x=element_blank()) +
# 	#	scale_x_discrete(labels=shortlab)
# 	# 	ylim(0,0.6)
# 	# print(g)
# 	# dev.off()

# }

# fst(1e-5)

# ## ASB ========================================================================================================================================

# output of gtex_v8_rsid_map.py, allele_specific_binding.py
asb_map <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/gtex_snp_to_asb_map.tsv'
asb_map <- fread(asb_map, header=FALSE, sep='\t')
colnames(asb_map) <- c('variant_id','asb_rsid')
asb_map <- asb_map[!is.na(asb_map[,asb_rsid])]

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

master_asb <- merge(egenes_master, asb_map, by='variant_id', all.x=TRUE)
master_asb[,annotation := ifelse(is.na(asb_rsid),0,1)]

# output of allele_specific_binding.py

asb <- function(cutoff){

	# dflist <- list()
	dfsame <- list()
	counts_df <- list()
	i <- 1
	for (t in tissues){

		dt <- master_asb[tissue==t]
		dt <- dt[pval_nominal < cutoff]

		# right now, this includes same SNPs
		dt <- dt[,.(gene_id, variant_id, pval_nominal, method, annotation)]
		dt[,pair := paste(gene_id, variant_id, sep=':')]

		global <- dt[method=='global']
		lava <- dt[method=='LAVA']

		# remove eGenes that have overlapping top SNPs
		repeats <- unique(lava[pair %in% global[,pair],gene_id])
		global <- global[!(gene_id %in% repeats)]
		lava <- lava[!(gene_id %in% repeats)]

		# get counts of genes that have >=1 annotated top SNP
		global <- global[,.(gene_id, annotation)]
		lava <- lava[,.(gene_id, annotation)]

		g_annot <- global[,list(summ_annot=sum(annotation)), by = gene_id]
		g_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
		g_annot[,summ_annot:=NULL]
		l_annot <- lava[,list(summ_annot=sum(annotation)), by = gene_id]
		l_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
		l_annot[,summ_annot:=NULL]

		# merge 
		all <- merge(g_annot, l_annot, by='gene_id', all=TRUE, suffixes=c('.global','.lava'))

		# # first, all eGenes where lead SNPs are different:
		# total_global <- nrow(g_annot)
		# total_lava <- nrow(l_annot)
		# g_hit <- sum(all[,bin.global], na.rm=TRUE)
		# l_hit <- sum(all[,bin.lava], na.rm=TRUE)

		# pval <- fisher(l_hit, g_hit, total_lava, total_global)

		# g_frac <- g_hit/total_global
		# l_frac <- l_hit/total_lava

		# df <- data.frame(tissue=t,
		# 	method=c('global','LAVA'),
		# 	frac=c(g_frac, l_frac),
		# 	count=c(total_global, total_lava),
		# 	pval=pval)
		# print(df)

		# dflist[[i]] <- df

		# now, only overlapping eGenes
		all <- all[complete.cases(all)]
		counts_df[[i]] <- all
		# tab <- table(all[,bin.global], all[,bin.lava])
		# if (nrow(tab)==2 & ncol(tab)==2){
		# 	a <- tab[1,1]
		# 	b <- tab[1,2]
		# 	c <- tab[2,1]
		# 	d <- tab[2,2]
		# 	pval <- mcnemar(a,b,c,d)
		# } else if (nrow(tab)==1 & ncol(tab)==2) {
		# 	a <- tab[1,1]
		# 	b <- tab[1,2]
		# 	pval <- mcnemar(a,b,0,0)
		# } else {
		# 	pval <- NA
		# }
		# print(paste(t, pval))

		total <- nrow(all)
		g_hit <- sum(all[,bin.global])
		l_hit <- sum(all[,bin.lava])

		df2 <- data.frame(tissue=t,
			method=c('global','LAVA'),
			hit=c(g_hit, l_hit),
			count=total)

		dfsame[[i]] <- df2

		i <- i + 1

	}

	# master1 <- data.table(rbindlist(dflist))
	# #print(master1)

	# label_pos <- max(master1[,frac]) + 0.05

	# g1 <- ggplot(master1, aes(x=factor(tissue),y=frac,fill=factor(method))) +
	# 	geom_bar(stat='identity', position='dodge') +
	# 	geom_text(aes(label = count, x = tissue, y = frac + 0.01), position = position_dodge(width = 0.9), vjust = -0.5) +
	# 	geom_text(aes(label = round(pval, 3), x = tissue, y = label_pos - 0.03), vjust = -0.5) +
	# 	scale_fill_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
	# 	labs(y='% ASB sites in top SNPs') +
	# 	theme(axis.title.x=element_blank()) +
	#	scale_x_discrete(labels=shortlab) +
	# 	ylim(0,label_pos)

	# pdf(paste0("~/gtex-admix/plots/final_figures/asb-all-",cutoff,".pdf"),width=10,height=4)
	# print(g1)
	# dev.off()

	master2 <- data.table(rbindlist(dfsame))
	labels <- master2[,.(tissue, count)]
	labels <- unique(labels)
	print(labels)

	label_pos <- 1.1*max(master2[,hit])

	g2 <- ggplot() +
		geom_bar(data=master2, aes(x=factor(tissue),y=hit,fill=factor(method)), stat='identity', position='dodge') +
		geom_text(data=labels, aes(label = paste0('N = ',count), x = tissue, y = label_pos), vjust = -0.5) +
		scale_fill_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='Local'), name='Method') +
		labs(y='N top SNP sets with ASB sites') +
		theme_bw() + 
		theme(axis.title.x=element_blank(),
			axis.text.x=element_text(colour='black',size=10)) +
		scale_x_discrete(labels=shortlab) +
		ylim(0,1.2*max(master2[,hit])) 

	pdf(paste0("~/gtex-admix/plots/final_figures/asb-same-",cutoff,".pdf"),width=10,height=4)
	print(g2)
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

# now do it for different cutoffs

# asb(1e-5)
# asb(1e-6)


# ## TF binding ========================================================================================================================================

# # output of encode_tf_binding.py

# tf_bind <- function(cutoff){

# 	#dflist <- list()
# 	dfsame <- list()
# 	counts_df <- list()
# 	i <- 1
# 	for (t in tissues){

# 		# define infiles 
# 		if (t == 'Muscle_Skeletal'){
# 			label <- 'gastrocnemius_medialis'
# 		} else if(t == 'Skin_Not_Sun_Exposed_Suprapubic'){
# 			label <- 'suprapubic_skin'
# 		} else if(t == 'Lung'){
# 			label <- 'lung'
# 		} else if(t == 'Nerve_Tibial'){
# 			label <- 'tibial_nerve'
# 		} else if(t == 'Artery_Tibial'){
# 			label <- 'tibial_artery'
# 		} else if(t == 'Adipose_Subcutaneous'){
# 			label <- 'subcutaneous_adipose_tissue'
# 		}

# 		infile <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/',label,'/',label,'.tf.binding.tsv')
# 		dt <- fread(infile, sep='\t', header=TRUE)
# 		dt <- dt[pval_nominal < cutoff]

# 		# right now, this includes same SNPs
# 		dt <- dt[,.(gene_id, variant_id, pval_nominal, method, annotation)]
# 		dt[,pair := paste(gene_id, variant_id, sep=':')]

# 		global <- dt[method=='global']
# 		lava <- dt[method=='LAVA']

# 		# remove eGenes that have overlapping top SNPs
# 		repeats <- unique(lava[pair %in% global[,pair],gene_id])
# 		global <- global[!(gene_id %in% repeats)]
# 		lava <- lava[!(gene_id %in% repeats)]

# 		# get counts of genes that have >= annotated top SNP
# 		global <- global[,.(gene_id, annotation)]
# 		lava <- lava[,.(gene_id, annotation)]

# 		g_annot <- global[,list(summ_annot=sum(annotation)), by = gene_id]
# 		g_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
# 		g_annot[,summ_annot:=NULL]
# 		l_annot <- lava[,list(summ_annot=sum(annotation)), by = gene_id]
# 		l_annot[, bin := ifelse(summ_annot > 0, 1, 0)]
# 		l_annot[,summ_annot:=NULL]

# 		# merge 
# 		all <- merge(g_annot, l_annot, by='gene_id', all=TRUE, suffixes=c('.global','.lava'))

# 		# # first, all eGenes where lead SNPs are different:
# 		# total_global <- nrow(g_annot)
# 		# total_lava <- nrow(l_annot)
# 		# g_hit <- sum(all[,bin.global], na.rm=TRUE)
# 		# l_hit <- sum(all[,bin.lava], na.rm=TRUE)

# 		# pval <- chi_sq(l_hit, g_hit, total_lava, total_global)

# 		# g_frac <- g_hit/total_global
# 		# l_frac <- l_hit/total_lava

# 		# df <- data.frame(tissue=t,
# 		# 	method=c('global','LAVA'),
# 		# 	hit=c(g_hit, l_hit),
# 		# 	count=c(total_global, total_lava),
# 		# 	pval=pval)

# 		# dflist[[i]] <- df

# 		# now, only overlapping eGenes
# 		all <- all[complete.cases(all)]
# 		counts_df[[i]] <- all
# 		tab <- table(all[,bin.global], all[,bin.lava])
# 		a <- tab[1,1]
# 		b <- tab[1,2]
# 		c <- tab[2,1]
# 		d <- tab[2,2]
# 		pval <- mcnemar(a,b,c,d)

# 		total <- nrow(all)
# 		g_hit <- sum(all[,bin.global])
# 		l_hit <- sum(all[,bin.lava])

# 		df2 <- data.frame(tissue=t,
# 			method=c('global','LAVA'),
# 			hit=c(g_hit, l_hit),
# 			count=total,
# 			pval=pval)

# 		dfsame[[i]] <- df2

# 		i <- i + 1

# 	}

# 	# master1 <- data.table(rbindlist(dflist))
# 	# print(master1)

# 	# label_pos <- max(master1[,frac]) + 0.05

# 	# g1 <- ggplot(master1, aes(x=factor(tissue),y=hit,fill=factor(method))) +
# 	# 	geom_bar(stat='identity', position='dodge') +
# 	# 	geom_text(aes(label = count, x = tissue, y = frac + 0.01), position = position_dodge(width = 0.9), vjust = -0.5) +
# 	# 	geom_text(aes(label = round(pval, 3), x = tissue, y = label_pos - 0.03), vjust = -0.5) +
# 	# 	scale_fill_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
# 	# 	labs(y='% TF binding site in top SNPs') +
# 	# 	theme(axis.title.x=element_blank()) +
# 	# 	ylim(0,label_pos)

# 	# pdf(paste0("~/gtex-admix/plots/final_figures/tf-binding-all-",cutoff,".pdf"),width=10,height=4)
# 	# print(g1)
# 	# dev.off()

# 	master2 <- data.table(rbindlist(dfsame))
# 	print(master2)

# 	label_pos <- 1.2*max(master2[,hit]) 

# 	g2 <- ggplot(master2, aes(x=factor(tissue),y=hit,fill=factor(method))) +
# 		geom_bar(stat='identity', position='dodge') +
# 		geom_text(aes(label = paste0('N = ',count), x = tissue, y = 1.1*max(master2[,hit])), vjust = -0.5) +
# 		scale_fill_manual(values=c(global='blue',LAVA='red'), labels=c(global='Global', LAVA='LAVA'), name='Method') +
# 		labs(y='N top SNP sets with TF binding sites') +
# 		theme_bw() +
# 		theme(axis.title.x=element_blank(),
#  			axis.text.x=element_text(colour='black',size=10)) +
# 		scale_x_discrete(labels=shortlab) +
# 		ylim(0,label_pos)

# 	pdf(paste0("~/gtex-admix/plots/final_figures/tf-binding-same-",cutoff,".pdf"),width=10,height=4)
# 	print(g2)
# 	dev.off()

# 	# get pval across all tissues 
# 	counts <- data.table(rbindlist(counts_df))
# 	counts <- unique(counts)

# 	tab <- table(counts[,bin.global], counts[,bin.lava])
# 	print(tab)
# 	a <- tab[1,1]
# 	b <- tab[1,2]
# 	c <- tab[2,1]
# 	d <- tab[2,2]
# 	pval <- mcnemar(a,b,c,d)
# 	print(pval)

# }

# # now do it for different cutoffs

# tf_bind(1e-6)


## coloc finemapped ========================================================================================================================================

# # snpmap <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/master/gtex_hg19/gtex_v7_to_v8_snps.txt'
# # finemap <- '/users/mgloud/projects/gwas/output/all_finemap_results_sorted_filtered.txt'
# # esnps <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/master/master.eSNP.txt.gz'

# # snps <- fread(cmd=sprintf('zcat %s', esnps), sep='\t', header=FALSE)
# # colnames(snps) <- 'variant_id'
# # split <- function(x){
# # 	return(paste(unname(unlist(strsplit(x,'_')))[1:2], collapse='_'))
# # }
# # snps[,v8 := sapply(variant_id, split)]

# # build_map <- fread(snpmap, sep='\t', header=TRUE)
# # snpmap <- merge(snps, build_map, by='v8')

# # get_tissue <- function(x){
# # 	return(unname(unlist(strsplit(x,'_allpairs')))[1])
# # }

# # get_trait <- function(x){
# # 	l = unname(unlist(strsplit(x,'/')))
# # 	return(l[length(l)])
# # }

# # finemap <- fread(finemap, sep='\t', header=TRUE)
# # finemap[,ref_snp := paste0('chr',ref_snp)]
# # finemap[,tissue := sapply(eqtl_file, get_tissue) ]
# # finemap[,phenotype := sapply(gwas_trait,get_trait) ]

# # master_fine <- merge(finemap, snpmap, by.x='ref_snp', by.y='v7')
# # master_fine <- master_fine[,.(variant_id, feature, tissue, phenotype)]
# # head(master_fine)

# # load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

# # master_fine <- merge(egenes_master, master_fine, by.x=c('gene_id','variant_id','tissue'), by.y=c('feature','variant_id','tissue'), all.x=TRUE)
# # head(master_fine)

# # save(master_fine, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/master_fine.v7.RData')

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/master_fine.v7.RData')

# master_fine[, hit := ifelse(is.na(phenotype), 0, 1)]
# table(master_fine[,hit])

# # only 16 pairs with SNPs that have hits in finemapping results, and none of them are marginally significant in my eQTL dataset 

# finemap <- function(cutoff){
# 	fine <- master_fine[pval_nominal < cutoff]

# 	# first, for all egenes with different lead snps
# 	fine[, pair := paste(gene_id,variant_id,sep=':')]
# 	lava <- fine[method=='LAVA']
# 	global <- fine[method=='global']
# 	repeats <- unique(lava[pair %in% global[,pair],gene_id])
# 	global <- global[!(gene_id %in% repeats)]
# 	lava <- lava[!(gene_id %in% repeats)]

# 	fine <- data.table(rbind(global, lava))
# 	print(table(fine[,hit]))

# 	hit <- as.data.frame.matrix(table(fine[hit==1,tissue],fine[hit==1,method]))
# 	hit <- data.table(hit, keep.rownames=TRUE)
# 	hit[,hit:=1]
# 	nohit <- as.data.frame.matrix(table(fine[hit==0,tissue],fine[hit==0,method]))
# 	nohit <- data.table(nohit, keep.rownames=TRUE)
# 	nohit[,hit:=0]

# 	counts <- data.table(rbind(hit,nohit))
# 	print(counts)
# }


## Biological example: T2D ========================================================================================================================================

# expected tissues: muscle, fat 

# rs231356:KCNQ1
# chr11:2684113 --> ENSG00000053918

# rs2244020:HLA-B
# chr6:31379674 --> ENSG00000234745

# rs231356
# rs2283228


# variant_id					tissue					gene_id				gene_name								pval_global				pval_lava
# chr11_2684113_A_T_b38		Adipose_Subcutaneous	ENSG00000238184.1	CD81 antisense RNA 1					0.0438148442538366		0.814449696269458
# chr11_2684113_A_T_b38		Muscle_Skeletal			NA					NA										NA						NA 
# chr11_2684113_A_T_b38		Artery_Tibial			ENSG00000229512.1	Novel transcript						0.0317674015382804		0.19518747879115
# chr11_2684113_A_T_b38		Artery_Tibial			ENSG00000236264.4	ribosomal protein L26 pseudogene 30		0.0313205501196163		0.912657940755304
# chr11_2684113_A_T_b38		Artery_Tibial			ENSG00000238184.1	CD81 antisense RNA 1					0.00220421167574171		0.875487969421866
# chr6_31379674				Adipose_Subcutaneous	NA					NA										NA						NA
# chr6_31379674				Muscle_Skeletal			NA					NA										NA						NA
# chr6_31379674				Artery_Tibial			NA					NA										NA						NA
# chr11_2828300_A_C_b38		Adipose_Subcutaneous	ENSG00000214026.10	mitochondrial ribosomal protein L23		0.042224073474703		0.513259969753511
# chr11_2828300_A_C_b38		Adipose_Subcutaneous	ENSG00000183734.4	achaete-scute family bHLH transcription factor 2	0.0469270576194101		0.901426985051925
# chr11_2828300_A_C_b38		Adipose_Subcutaneous	ENSG00000184281.14	tumor suppressing subtransferable candidate 4	0.0365162097117613		0.244078667316157
# chr11_2828300_A_C_b38		Adipose_Subcutaneous	ENSG00000021762.19	oxysterol binding protein like 5		0.0283652381362871		0.771337687244138
# chr11_2828300_A_C_b38		Adipose_Subcutaneous	ENSG00000177105.9	ras homolog family member G				0.0555551911277125		0.130986426487654


