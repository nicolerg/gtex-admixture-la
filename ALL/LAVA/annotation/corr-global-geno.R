library(data.table)
library(ggplot2)

# See version with more SNPs: r2_all_snps.R

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

# show correlation between genotypes and local ancestry at global eQTL loci

# first, get a list of global top SNPs (one from each eGene)
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')
head(egenes_master)
cutoff <- 1e-5
egenes_master <- egenes_master[pval_nominal < cutoff]
# get matching egenes

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

	# remove eGenes that are not in both methods
	global <- global[gene_id %in% lava[,gene_id]]
	lava <- lava[gene_id %in% global[,gene_id]]

	# # pick first variant for each gene
	# global <- global[!duplicated(gene_id)]
	# lava <- lava[!duplicated(gene_id)]

	# need to keep track of how many values there are for each tissue and method 
	m <- data.table(rbind(global,lava))

	dtlist[[i]] <- m
	i <- i+1
}

master <- data.table(rbindlist(dtlist))
print(table(master[,tissue],master[,method]))

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
# save(all_snps, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cor/all_snps.RData')

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cor/all_snps.RData')

mas <- merge(master, all_snps, by='variant_id')

mas[gene_id == "ENSG00000224950.2" & tissue == 'Adipose_Subcutaneous']
head(mas)
q()

dtlist <- list()
countdf <- list()
i<-1 
for (t in tissues){

	dt <- mas[tissue==t]

	global <- dt[method=='global']
	lava <- dt[method=='LAVA']

	# remove eGenes that have overlapping top SNPs
	repeats <- unique(lava[pair %in% global[,pair],gene_id])
	global <- global[!(gene_id %in% repeats)]
	lava <- lava[!(gene_id %in% repeats)]

	# remove eGenes that are not in both methods
	global <- global[gene_id %in% lava[,gene_id]]
	lava <- lava[gene_id %in% global[,gene_id]]

	# keep only those whose difference in p-value is above some cutoff 

	global <- global[,.(gene_id, pval_nominal, tissue, local_rsq)]
	lava <- lava[,.(gene_id, pval_nominal, tissue, local_rsq)]

	global <- global[,list(max_score=max(local_rsq)), by = .(gene_id,tissue,pval_nominal)]
	lava <- lava[,list(max_score=max(local_rsq)), by = .(gene_id,tissue,pval_nominal)]

	m <- merge(global, lava, by=c('gene_id','tissue'), suffixes=c('.global','.lava'))
	m_filt <- m[abs(-log10(pval_nominal.global) - -log10(pval_nominal.lava))>=1]

	dtlist[[i]] <- m_filt
	i <- i+1

	print(t.test(m_filt[,max_score.lava], m_filt[,max_score.global]))

}

master_cor <- data.table(rbindlist(dtlist))
counts <- as.data.table(table(master_cor[,tissue]))
colnames(counts) <- c('tissue', 'N')

print(t.test(master_cor[,max_score.lava], master_cor[,max_score.global]))

master_cor[max_score.global==max(max_score.global)]

master_cor <- melt(master_cor, measure.vars=c('max_score.global', 'max_score.lava'))


pdf(paste0("~/gtex-admix/plots/final_figures/local-rsq-",cutoff,"-same.pdf"),width=8,height=4)
g<-ggplot() +
	geom_boxplot(data=master_cor,aes(x=factor(tissue),y=value,fill=factor(variable))) +
	geom_text(data=counts, aes(x=factor(tissue), y=1, label=paste0('N = ',N)), colour='black') +
	scale_fill_manual(values = c('max_score.global'='blue', 'max_score.lava'='red'),
		labels=c('max_score.lava'='Local','max_score.global'='Global'), name='Method') +
	labs(y='R-squared between genotype and local ancestry') +
	theme_bw() +
	theme(axis.title.x=element_blank(),
		axis.text.x=element_text(colour='black')) +
	scale_x_discrete(labels=shortlab)
print(g)
dev.off()

