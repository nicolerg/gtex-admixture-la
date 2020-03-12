library(data.table)
library(ggplot2)
library(ggrepel)

shortlab <- c(Adipose_Subcutaneous='Adipose',
	Artery_Tibial='Tibial\nartery',
	Lung='Lung',
	Muscle_Skeletal='Skeletal\nmuscle',
	Nerve_Tibial='Tibial\nnerve',
	Skin_Not_Sun_Exposed_Suprapubic='Skin (not\nsun-exposed)')

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

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData') # top SNPs per gene 

cutoff <- 1e-6
egenes_master <- egenes_master[pval_nominal < cutoff]

dtlist <- list()
i<-1 
for (t in tissues){

	dt <- egenes_master[tissue==t]
	dt[,pair := paste(gene_id, variant_id, sep=':')]

	global <- dt[method=='global']
	lava <- dt[method=='LAVA']

	# # remove eGenes that have overlapping top SNPs
	# repeats <- unique(lava[pair %in% global[,pair],gene_id])
	# global <- global[!(gene_id %in% repeats)]
	# lava <- lava[!(gene_id %in% repeats)]

	# unique eGenes only 
	g <- global[!(gene_id %in% lava[,gene_id])]
	l <- lava[!(gene_id %in% global[,gene_id])]

	# # pick first variant for each gene
	# global <- global[!duplicated(gene_id)]
	# lava <- lava[!duplicated(gene_id)]

	m <- data.table(rbind(g,l))

	dtlist[[i]] <- m
	i <- i+1
}

master <- data.table(rbindlist(dtlist))


# DEPRECATED - see r2_all_snps.R ===================================================
# get_chrom <- function(x){
# 	return(unname(unlist(strsplit(x, '_')))[1])
# }

# snps_of_interest <- unique(master[,variant_id])
# length(snps_of_interest)
# snps <- data.table(variant_id=snps_of_interest,
# 	chromosome=sapply(snps_of_interest, get_chrom))
# snps <- snps[order(chromosome)]

# # get local ancestry for each of those SNPs 
# cov_base='/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/'
# vcf_base='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/'
# outdir='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cor/'

# snp_dt <- list()

# get_n_hap <- function(x){
# 	n <- unname(unlist(strsplit(x, "\\|")))
# 	if (n[1] == '.'){
# 		return(NA)
# 	} else {
# 		return(as.integer(n[1])+as.integer(n[2]))
# 	}
# }

# convert_col <- function(x){
# 	return(sapply(x, get_n_hap))
# }

# i <- 1
# for (chr in unique(snps[,chromosome])){
# 	snp_list <- snps[chromosome == chr, variant_id]

# 	# get local ancestry 
# 	local_cov <- paste0(cov_base, chr, '.hg19.maf0.localcov.tsv.gz')
# 	system('echo "Reading in local cov file"')
# 	local <- fread(cmd=sprintf("zcat %s", local_cov), header=TRUE, sep='\t')
# 	# SUBJID  SNP_ID  POS     AFR     ASN     EUR     UNK
# 	local <- local[SNP_ID %in% snp_list, .(SUBJID, SNP_ID, AFR, ASN, EUR)]
# 	system('echo "Filtering local cov file"')
# 	print(head(local))

# 	subjs <- unique(local[,SUBJID])

# 	# get genotypes for each SNP
# 	vcf <- paste0(vcf_base, 'gtex.admixed.MAC10.phased.', chr, '.recode.vcf.gz')
# 	out <- paste0(outdir, chr, '.snps')
# 	cmd <- sprintf('zcat %s | head -1 > %s', vcf, out)
# 	system(cmd)
# 	snp_file <- data.table(id=snp_list)
# 	snp_path <- paste0(outdir,chr,'.snps.txt')
# 	write.table(snp_file, snp_path, sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
# 	system('echo "Generating filtered VCF"')
# 	cmd <- sprintf('vcftools --gzvcf %s --out %s --snps %s --recode', vcf, out, snp_path)
# 	system(cmd)

# 	system('echo "Transforming VCF"')
# 	vcf_filt <- fread(paste0(out,".recode.vcf"), sep='\t', header=TRUE)
# 	cols <- subjs[subjs %in% colnames(vcf_filt)]
# 	vcf_filt <- vcf_filt[,c('ID', cols),with=FALSE]
# 	vcf_filt[ , (cols) := lapply(.SD, convert_col), .SDcols = cols]

# 	row_label <- vcf_filt[,ID]
# 	vcf_filt[,ID := NULL]
# 	subjid <- colnames(vcf_filt)
# 	vcf_filt <- data.table(t(vcf_filt))
# 	colnames(vcf_filt) <- row_label
# 	vcf_filt[,SUBJID := subjid]

# 	c <- data.table(variant_id=snp_list, local_rsq=0.0)
# 	for (v in snp_list){
# 		g <- vcf_filt[,c('SUBJID', v), with=FALSE]
# 		l <- local[SNP_ID == v]
# 		m <- merge(g, l, by='SUBJID')
# 		f <- as.formula(sprintf('%s ~ AFR + ASN',v))
# 		lm.fit <- lm(f, data=m)
# 		rsq <- summary(lm.fit)$r.squared
# 		c[variant_id==v, local_rsq := rsq]
# 	}

# 	snp_dt[[i]] <- c
# 	i <- i + 1 

# }

# all_snps <- data.table(rbindlist(snp_dt))
# nrow(all_snps)
# head(all_snps)
# save(all_snps, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cor/all_snps-diff.RData')
# DEPRECATED - see r2_all_snps.R ===================================================

# ==============================================================================================================================
# # unique eGenes only (at given cutoff)

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cor/all_snps-diff.RData')
# mas <- merge(master, all_snps, by='variant_id', all.x=TRUE)

# dtlist <- list()
# countdf <- list()
# i<-1 
# for (t in tissues){

# 	dt <- mas[tissue==t]

# 	global <- dt[method=='global']
# 	lava <- dt[method=='LAVA']

# 	# # # remove eGenes that have overlapping top SNPs
# 	# # repeats <- unique(lava[pair %in% global[,pair],gene_id])
# 	# # global <- global[!(gene_id %in% repeats)]
# 	# # lava <- lava[!(gene_id %in% repeats)]

# 	# remove eGenes that are in both methods
# 	global <- global[!(gene_id %in% lava[,gene_id])]
# 	lava <- lava[!(gene_id %in% global[,gene_id])]

# 	global <- global[,.(gene_id, pval_nominal, tissue, local_rsq)]
# 	lava <- lava[,.(gene_id, pval_nominal, tissue, local_rsq)]

# 	global <- global[,list(max_score=max(local_rsq)), by = .(gene_id,tissue,pval_nominal)]
# 	lava <- lava[,list(max_score=max(local_rsq)), by = .(gene_id,tissue,pval_nominal)]
# 	global[,method:='global']
# 	lava[,method:='LAVA']

# 	m <- data.table(rbind(global, lava))

# 	dtlist[[i]] <- m
# 	i <- i+1

# }

# master_cor <- data.table(rbindlist(dtlist))

# counts <- as.data.table(table(master_cor[,tissue], master_cor[,method]))
# colnames(counts) <- c('tissue', 'method', 'N')
# counts

# pdf(paste0("~/gtex-admix/plots/final_figures/local-rsq-",cutoff,"-diff.pdf"),width=8,height=4)
# g<-ggplot() +
# 	geom_boxplot(data=master_cor,aes(x=factor(tissue),y=max_score,fill=factor(method))) +
# 	geom_text(data=counts, aes(x=factor(tissue), y=1, label=N, colour=factor(method)), position=position_dodge(width=0.6)) +
# 	scale_colour_manual(values = c('global'='blue', 'LAVA'='red'),
# 		labels=c('LAVA'='Local','global'='Global'), name='Method') +
# 	scale_fill_manual(values = c('global'='blue', 'LAVA'='red'),
# 		labels=c('LAVA'='Local','global'='Global'), name='Method') +
# 	labs(y='R-squared between genotype and local ancestry') +
# 	theme_bw() +
# 	theme(axis.title.x=element_blank(),
# 		axis.text.x=element_text(colour='black')) +
# 	scale_x_discrete(labels=shortlab)
# print(g)
# dev.off()

# ==============================================================================================================================

# 1. look at p-values for identical tests in local (consider unique global eSNPs only)

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cor/all_snps-diff.RData')
mas <- merge(master, all_snps, by='variant_id', all.x=TRUE)
gene_labels <- fread('~/gtex-admix/metadata/gene_id_conversion.txt',header=TRUE,sep='\t')
head(gene_labels)

# plot p-value difference for unique eQTLs in each method
for (m in c('global','LAVA')){
	dt <- mas[method==m]

	outfile=paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/dt-diff',m,'-match.RData')
	if(file.exists(outfile)){
		load(outfile)
	} else {
		if(!exists('allp_master')){
			load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/allp_master-05.RData')
			allp_master <- allp_master[,.(gene_id, variant_id, tissue, global_pval_nominal, global_slope, lava_pval_nominal, lava_slope)]
		}
		dt <- merge(dt, allp_master, by=c('gene_id', 'variant_id', 'tissue'), all.x=TRUE)
		dt[,method:=NULL]
		dt[,pair:=NULL]
		dt[,pval_nominal:=NULL]
		dt[,slope:=NULL]	
		save(dt,file=outfile)
	}

	if(m=='global'){
		dt[,diff := -log10(global_pval_nominal) + log10(lava_pval_nominal)]
	} else {
		dt[,diff := -log10(lava_pval_nominal) + log10(global_pval_nominal)]
	}

	dt <- dt[order(diff, decreasing=TRUE)]
	dt <- merge(dt, gene_labels, by.x='gene_id', by.y='Gene stable ID version',all.x=TRUE)

	plot_df <- dt[,.(global_pval_nominal, lava_pval_nominal, diff, gene_id, `Gene name`, tissue)]
	plot_df <- plot_df[!duplicated(plot_df)]
	plot_df[,label := ifelse(diff < 3, NA, `Gene name`)]

	if(m=='global'){
		plot_df[,x_start := 'Global']
		plot_df[,x_end := 'Local']

		pdf("~/gtex-admix/plots/final_figures/unique-global-into-local-pval.pdf",width=10,height=8)
		g<-ggplot(plot_df, aes(colour=factor(tissue))) +
			geom_segment(aes(x=x_start, xend=x_end, y=-log10(global_pval_nominal), yend=-log10(lava_pval_nominal)), alpha=0.5) +
			geom_text_repel(aes(x=x_end, y=-log10(lava_pval_nominal), label=label), hjust=0, direction='y', nudge_x=1, size=3) +
			labs(y='-log10(pvalue)',title='Unique Global eQTLs into Local tests') +
			theme_bw() +
			theme(axis.title.x=element_blank()) +
			scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue')
		print(g)
		dev.off()

	} else {
		plot_df[,x_start := 'Local']
		plot_df[,x_end := 'Global']

		pdf("~/gtex-admix/plots/final_figures/unique-local-into-global-pval.pdf",width=10,height=8)
		g<-ggplot(plot_df, aes(colour=factor(tissue))) +
			geom_segment(aes(x=x_start, xend=x_end, y=-log10(lava_pval_nominal), yend=-log10(global_pval_nominal)), alpha=0.5) +
			geom_text_repel(aes(x=x_end, y=-log10(global_pval_nominal), label=label), hjust=0, direction='y', nudge_x=1, size=3) +
			labs(y='-log10(pvalue)', title='Unique Local eQTLs into Global tests') +
			theme_bw() +
			scale_x_discrete(breaks=c('Local','Global'),limits=c('Local','Global')) +
			theme(axis.title.x=element_blank()) +
			scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue')
		print(g)
		dev.off()

	}

}



# compare to genes in Mogil

pop_gene_r2 <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/Mogil/r2/merged_r2_mogil.tsv'























# # plot p-value different by R^2 between genotype and local ancestry at unique eSNPs
# for (m in c('global','LAVA')){
# 	dt <- mas[method==m]

# 	outfile=paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/dt-diff',m,'-match.RData')
# 	if(file.exists(outfile)){
# 		load(outfile)
# 	} else {
# 		if(!exists('allp_master')){
# 			load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/allp_master-05.RData')
# 			allp_master <- allp_master[,.(gene_id, variant_id, tissue, global_pval_nominal, global_slope, lava_pval_nominal, lava_slope)]
# 		}
# 		dt <- merge(dt, allp_master, by=c('gene_id', 'variant_id', 'tissue'), all.x=TRUE)
# 		dt[,method:=NULL]
# 		dt[,pair:=NULL]
# 		dt[,pval_nominal:=NULL]
# 		dt[,slope:=NULL]	
# 		save(dt,file=outfile)
# 	}

# 	if(m=='global'){
# 		dt[,diff := -log10(global_pval_nominal) + log10(lava_pval_nominal)]
# 	} else {
# 		dt[,diff := -log10(lava_pval_nominal) + log10(global_pval_nominal)]
# 	}

# 	dt <- dt[order(diff, decreasing=TRUE)]

# 	plot_df <- dt[,.(local_rsq, diff)]

# 	if(m=='global'){

# 		pdf("~/gtex-admix/plots/final_figures/unique-global-diff-vs-rsq.pdf",width=4,height=4)
# 		g<-ggplot(plot_df, aes(x=diff, y=local_rsq)) +
# 			geom_point() +
# 			labs(x='-log10(Global) - -log10(Local)', y='R-squared between genotype and local ancestry') +
# 			theme_bw() 
# 		print(g)
# 		dev.off()

# 	} else {

# 		pdf("~/gtex-admix/plots/final_figures/unique-local-diff-vs-rsq.pdf",width=4,height=4)
# 		g<-ggplot(plot_df, aes(x=diff, y=local_rsq)) +
# 			geom_point() +
# 			labs(x='-log10(Local) - -log10(Global)', y='R-squared between genotype and local ancestry') +
# 			theme_bw() 
# 		print(g)
# 		dev.off()
# 	}

# }

# 2. look at R-squared for corresponding lead SNP in local (consider unique global eSNPs only)

# get unique global at cutoff of 1e-6
# go into local and get corresponding top SNPs for those egenes
# find correlation with local ancestry
# make line segment plot 

# I can do this once the other script is done running 






