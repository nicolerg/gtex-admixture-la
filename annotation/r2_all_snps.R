library(data.table)
library(ggplot2)
library(OneR)
library(plyr)

# calculate variance in genotype explained by local ancestry for my lead SNPs 
# use only 117AX genotypes and local ancestry 

args <- commandArgs(trailingOnly=TRUE)
chrom <- args[1]

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic',
	'Whole_Blood')

outdir='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/REVISIONS/'
cov_base='/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/'
vcf_base='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/'

snp_path=paste0(outdir,'/snp_list.txt')
if (!file.exists(snp_path)){

	# first, get SNP subset 
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/merged/egenes_master.RData')
	variants <- unique(egenes_master[tissue%in%tissues,variant_id])
	
	# for file in `ls /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL | grep "egenes"`; do zcat $file | cut -f12 | sed -e '1d' >> /mnt/lab_data/montgomery/nicolerg/gtex_esnp_list.txt; done
	# gtex_esnps <- fread(paste0(outdir,'/gtex.v8.esnp.txt'),sep='\t',header=FALSE)
	# variants <- c(variants, gtex_esnps[,V1])
	
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

	snps = fread(snp_path, sep='\t', header=F)
	snp_list = snps[,V1]

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

get_r2(chrom)
