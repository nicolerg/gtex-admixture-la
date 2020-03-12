library(data.table)
library(ggplot2)
library(OneR)
library(plyr)

# stratify R^2 by difference in EUR and AFR allele frequency 
# all SNPs in egenes_master; all lead SNPs from GTEx v8 

args <- commandArgs(trailingOnly=TRUE)
chrom <- args[1]
if (length(args) > 1){
	run_get_r2 <- as.numeric(args[2])
} else {
	run_get_r2 <- 0
}


tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

outdir='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/'
cov_base='/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/'
vcf_base='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/'

snp_path=paste0(outdir,'/snp_list.txt')
if (!file.exists(snp_path)){

	# first, get SNP subset 
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')
	variants <- unique(egenes_master[,variant_id])
	# for file in `ls /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL | grep "egenes"`; do zcat $file | cut -f12 | sed -e '1d' >> /mnt/lab_data/montgomery/nicolerg/gtex_esnp_list.txt; done
	gtex_esnps <- fread(paste0(outdir,'/gtex.v8.esnp.txt'),sep='\t',header=FALSE)
	variants <- c(variants, gtex_esnps[,V1])
	snp_list <- unique(variants)

	snp_df <- data.frame(variant_id=snp_list)
	write.table(snp_df,snp_path,sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
}

get_n_hap <- function(x){
	n <- unname(unlist(strsplit(x, "\\|")))
	if (n[1] == '.'){
		return(NA)
	} else {
		return(as.integer(n[1])+as.integer(n[2]))
	}
}

convert_col <- function(x){
	return(sapply(x, get_n_hap))
}

get_r2 <- function(chrom){

	# read in local ancestry
	local <- paste0(cov_base,'/chr',chrom,'.hg19.maf0.localcov.tsv.gz')
	local <- fread(cmd=sprintf('zcat %s',local), sep='\t', header=TRUE)
	local <- local[,.(SUBJID, SNP_ID, AFR, ASN, EUR)]
	subjs <- unique(local[,SUBJID])
	# filter SNPs
	local <- local[SNP_ID %in% snp_list]

	# get genotypes for each SNP
	vcf <- paste0(vcf_base, 'gtex.admixed.MAC10.phased.chr', chrom, '.recode.vcf.gz')
	out <- paste0(outdir, chrom, '.snps')
	if (!file.exists(paste0(out,".recode.vcf"))){
		# cmd <- sprintf('zcat %s | head -1 > %s', vcf, out)
		# system(cmd)
		system('echo "Generating filtered VCF"')
		cmd <- sprintf('vcftools --gzvcf %s --out %s --snps %s --recode', vcf, out, snp_path)
		system(cmd)
	}

	system('echo "Transforming VCF"')
	vcf_filt <- fread(paste0(out,".recode.vcf"), sep='\t', header=TRUE)
	cols <- subjs[subjs %in% colnames(vcf_filt)]
	vcf_filt <- vcf_filt[,c('ID', cols),with=FALSE]
	vcf_filt[ , (cols) := lapply(.SD, convert_col), .SDcols = cols]

	row_label <- vcf_filt[,ID]
	vcf_filt[,ID := NULL]
	subjid <- colnames(vcf_filt)
	vcf_filt <- data.table(t(vcf_filt))
	colnames(vcf_filt) <- row_label
	vcf_filt[,SUBJID := subjid]

	# create a progress bar
	total <- length(snp_list)
	count = 0
	pb <- txtProgressBar(min = 0, max = total, style = 3)

	outfile <- file(paste0(outdir,chrom,'.r2.txt'), 'w')
	cat("variant_id\tlocal_rsq\n", file=outfile)

	for (v in snp_list){
		
		# skip SNPs not in this file
		if (!v %in% colnames(vcf_filt)){
			next
		}

		g <- vcf_filt[,c('SUBJID', v), with=FALSE]
		l <- local[SNP_ID == v]
		m <- merge(g, l, by='SUBJID')
		f <- as.formula(sprintf('%s ~ AFR + ASN',v))
		lm.fit <- lm(f, data=m)
		rsq <- summary(lm.fit)$r.squared

		string <- paste(v, rsq, sep='\t')
		cat(paste0(string,'\n'), file=outfile, append=TRUE)

		count <- count + 1
		setTxtProgressBar(pb, count)
	}
	close(pb)
	close(outfile)
}

# get R^2 for random non-lead SNPs with a small difference in EUR and AFR allele frequency 

get_r2_rand <- function(chrom){

	snp_path <- paste0(outdir,'/gtex.admixed.MAC10.random.52133.af.lt.0.18.snps')

	snps <- fread(snp_path,header=FALSE,sep='\t')
	colnames(snps) <- 'variant_id'
	snp_list <- snps[,variant_id]

	# read in local ancestry
	local <- paste0(cov_base,'/chr',chrom,'.hg19.maf0.localcov.tsv.gz')
	local <- fread(cmd=sprintf('zcat %s',local), sep='\t', header=TRUE)
	local <- local[,.(SUBJID, SNP_ID, AFR, ASN, EUR)]
	subjs <- unique(local[,SUBJID])
	# filter SNPs
	local <- local[SNP_ID %in% snp_list]

	# get genotypes for each SNP
	vcf <- paste0(vcf_base, 'gtex.admixed.MAC10.phased.chr', chrom, '.recode.vcf.gz')
	out <- paste0(outdir, chrom, '.snps')
	if (!file.exists(paste0(out,".recode.vcf"))){
		# cmd <- sprintf('zcat %s | head -1 > %s', vcf, out)
		# system(cmd)
		system('echo "Generating filtered VCF"')
		cmd <- sprintf('vcftools --gzvcf %s --out %s --snps %s --recode', vcf, out, snp_path)
		system(cmd)
	}

	system('echo "Transforming VCF"')
	vcf_filt <- fread(paste0(out,".recode.vcf"), sep='\t', header=TRUE)
	cols <- subjs[subjs %in% colnames(vcf_filt)]
	vcf_filt <- vcf_filt[,c('ID', cols),with=FALSE]
	vcf_filt[ , (cols) := lapply(.SD, convert_col), .SDcols = cols]

	row_label <- vcf_filt[,ID]
	vcf_filt[,ID := NULL]
	subjid <- colnames(vcf_filt)
	vcf_filt <- data.table(t(vcf_filt))
	colnames(vcf_filt) <- row_label
	vcf_filt[,SUBJID := subjid]

	# create a progress bar
	total <- length(snp_list)
	count = 0
	pb <- txtProgressBar(min = 0, max = total, style = 3)

	outfile <- file(paste0(outdir,chrom,'.rand.r2.txt'), 'w')
	cat("variant_id\tlocal_rsq\n", file=outfile)

	for (v in snp_list){
		
		# skip SNPs not in this file
		if (!v %in% colnames(vcf_filt)){
			next
		}

		g <- vcf_filt[,c('SUBJID', v), with=FALSE]
		l <- local[SNP_ID == v]
		m <- merge(g, l, by='SUBJID')
		f <- as.formula(sprintf('%s ~ AFR + ASN',v))
		lm.fit <- lm(f, data=m)
		rsq <- summary(lm.fit)$r.squared

		string <- paste(v, rsq, sep='\t')
		cat(paste0(string,'\n'), file=outfile, append=TRUE)

		count <- count + 1
		setTxtProgressBar(pb, count)
	}
	close(pb)
	close(outfile)
}


merge_r2 <- function(cutoff=1e-4){

	# merge non-random SNPs. split between my lead SNPs and GTEx v8 lead SNPs

	if(file.exists(paste0(outdir,'/m-merged.RData'))){
		load(paste0(outdir,'/m-merged.RData'))
	} else {

		if(file.exists(paste0(outdir,'/r2_af.RData'))){
			load(paste0(outdir,'/r2_af.RData'))
		} else {
			af = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/af/gtex.admixed.MAC10.snps.af.txt.gz'
			af <- fread(cmd=sprintf('zcat %s', af), sep='\t', header=TRUE)
			r2 = paste0(outdir,'/all_lead_snp_r2.txt.gz')
			r2 <- fread(cmd=sprintf('zcat %s', r2), sep='\t', header=TRUE)

			r2_af <- merge(af, r2, by='variant_id', all.y=TRUE)
			r2_af <- r2_af[,.(variant_id,EAS_AF,EUR_AF,AFR_AF,local_rsq)]
			r2_af <- unique(r2_af, by='variant_id')

			save(r2_af, file=paste0(outdir,'/r2_af.RData'))	
		}

		# first, for my lead SNPs 

		load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_',cutoff,'.RData')) # top SNPs per gene 

		m <- merge(diff_lead, r2_af, by=c('variant_id'), all.x=TRUE)
		m[,pair := NULL]
		m[,slope := NULL]
		m[,count := NULL]
		m[,local_rsq := as.numeric(local_rsq)]

		subset <- m[,list(local_rsq=max(local_rsq, na.rm=TRUE)), by = .(gene_id,tissue,method)]
		m <- merge(m, subset, by=c('gene_id','tissue','method','local_rsq'))
		m <- unique(m, by=c('gene_id','tissue','method','local_rsq'))
		print(head(m))

		# one line per gene
		global <- m[method=='global']
		local <- m[method=='LAVA']
		global[,method := NULL]
		local[,method := NULL]
		m <- merge(global, local, by=c('gene_id','tissue'), suffixes=c('.global','.local'))
		m[,eur_afr_af_diff := abs(as.numeric(EUR_AF.global) - as.numeric(AFR_AF.global))]
		m[,eur_asn_af_diff := abs(as.numeric(EUR_AF.global) - as.numeric(EAS_AF.global))]
		m[,afr_asn_af_diff := abs(as.numeric(AFR_AF.global) - as.numeric(EAS_AF.global))]

		save(m, file=paste0(outdir,'/m-merged.RData'))

		# next, for GTEx lead SNPs (I have to get all of those allele frequencies first)
	}

	# merge files for random SNPs

	af <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/af/gtex.admixed.MAC10.random.52133.af.lt.0.18.txt'
	r2 <- paste0(outdir,'/r2_rand.txt')

	af <- fread(af, sep='\t', header=TRUE)
	r2 <- fread(r2, sep='\t', header=TRUE)

	r2_af_rand <- merge(af, r2, by='variant_id', all.y=TRUE)
	r2_af <- r2_af_rand[,.(variant_id,EAS_AF,EUR_AF,AFR_AF,local_rsq)]
	r2_af <- unique(r2_af, by='variant_id')

	# remove lead SNPs 
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData') # top SNPs per gene 
	r2_af <- r2_af[!variant_id %in% unique(egenes_master[,variant_id])]

	# merge with previous results

	colnames(r2_af) <- paste0(colnames(r2_af),'.rand')
	r2_af[,eur_afr_af_diff := abs(EUR_AF.rand - AFR_AF.rand)]
	r2_af[,eur_afr_af_diff := abs(EUR_AF.rand - AFR_AF.rand)]

	r2_af[,eur_afr_af_diff := abs(as.numeric(EUR_AF.rand) - as.numeric(AFR_AF.rand))]
	r2_af[,eur_asn_af_diff := abs(as.numeric(EUR_AF.rand) - as.numeric(EAS_AF.rand))]
	r2_af[,afr_asn_af_diff := abs(as.numeric(AFR_AF.rand) - as.numeric(EAS_AF.rand))]

	master <- data.table(rbind.fill(r2_af, m))

	save(master, file=paste0(outdir,'/master-admixed-all-merged.RData'))

}

if (run_get_r2 == 1){
	get_r2(chrom)
	#get_r2_rand(chrom)
} else {
	merge_r2()
}





# # plot difference in allele frequency by R^2

# # af = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/af/gtex.admixed.MAC10.snps.af.txt.gz'
# # af <- fread(cmd=sprintf('zcat %s', af), sep='\t', header=TRUE)
# # r2 = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/allsnps.01.r2.txt.gz'
# # r2 <- fread(cmd=sprintf('zcat %s', r2), sep='\t', header=TRUE)

# # r2_af <- merge(af, r2, by='variant_id', all.y=TRUE)
# # r2_af <- r2_af[,.(variant_id,EAS_AF,EUR_AF,AFR_AF,local_rsq)]
# # r2_af <- unique(r2_af, by='variant_id')

# # save(r2_af, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/r2_af.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/r2_af.RData')
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData') # top SNPs per gene 

# cutoff <- 0.01
# egenes_master <- egenes_master[pval_nominal < cutoff]

# dtlist <- list()
# i<-1 
# for (t in tissues){

# 	dt <- egenes_master[tissue==t]
# 	dt[,pair := paste(gene_id, variant_id, sep=':')]

# 	global <- dt[method=='global']
# 	lava <- dt[method=='LAVA']

# 	# remove eGenes that have overlapping top SNPs
# 	repeats <- unique(lava[pair %in% global[,pair],gene_id])
# 	global <- global[!(gene_id %in% repeats)]
# 	lava <- lava[!(gene_id %in% repeats)]

# 	# remove eGenes that are not in both methods
# 	g <- global[gene_id %in% lava[,gene_id]]
# 	l <- lava[gene_id %in% global[,gene_id]]

# 	m <- data.table(rbind(g,l))

# 	dtlist[[i]] <- m
# 	i <- i+1
# }

# master <- data.table(rbindlist(dtlist))
# nrow(master)

# m <- merge(master, r2_af, by=c('variant_id'), all.x=TRUE)
# m[,pair := NULL]
# m[,slope := NULL]
# m[,gtex_cutoff := NULL]
# m[,count := NULL]
# m[,local_rsq := as.numeric(local_rsq)]

# subset <- m[,list(local_rsq=max(local_rsq, na.rm=TRUE)), by = .(gene_id,tissue,method)]
# m <- merge(m, subset, by=c('gene_id','tissue','method','local_rsq'))
# m <- unique(m, by=c('gene_id','tissue','method','local_rsq'))
# m[,local_corr := sqrt(as.numeric(local_rsq))]
# head(m)

# # one line per gene
# global <- m[method=='global']
# local <- m[method=='LAVA']
# global[,method := NULL]
# local[,method := NULL]
# m <- merge(global, local, by=c('gene_id','tissue'), suffixes=c('.global','.local'))
# m[,eur_afr_af_diff := abs(as.numeric(EUR_AF.global) - as.numeric(AFR_AF.global))]

# # save(m, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/m-merged.RData')

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/m-merged.RData')

# # bin difference in allele frequency
# m <- m[complete.cases(m)]
# m[,bin := bin(eur_afr_af_diff, nbins=6, method='length')]
# bin_n <- data.table(table(m[,bin]))
# colnames(bin_n) <- c('bin','label')

# for (b in bin_n[,bin]){
	
# 	bin_n[bin == b, pvalue := signif(t.test(m[bin == b, local_corr.local], m[bin == b, local_corr.global])$p.value, digits=3) ]

# }

# m <- melt(m, measure=c('local_corr.global','local_corr.local'))

# cutoff <- 0.01

# jpeg(paste0('~/gtex-admix/plots/final_figures/la_corr_stratified_by_delta_af_',cutoff,'.jpg'),res=150,height=800,width=1600)
# ggplot() +
# 	geom_boxplot(data=m, aes(x=bin, y=value, colour=variable)) +
# 	scale_colour_manual(values=c(local_corr.global='blue',local_corr.local='red'),
# 		labels=c(local_corr.global='Global',local_corr.local='Local'), name='Method') +
# 	labs(x='| Global lead SNP African AF - Global lead SNP European AF |',y='Correlation between lead SNP genotypes and local ancestry') +
# 	theme_bw() +
# 	geom_label(data=bin_n, aes(x=bin, y=max(m[,value]+0.12), label=paste0('N = ',label))) +
# 	geom_label(data=bin_n, aes(x=bin, y=max(m[,value]+0.05), label=pvalue))
# dev.off()
