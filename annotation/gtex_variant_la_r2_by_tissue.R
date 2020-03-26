# 9/18/2019

# find VE between GTEx v8 lead SNP genotype and LA, accounting for individuals in each tissue

library(data.table)
library(ggplot2)
library(OneR)
library(plyr)

args <- commandArgs(trailingOnly=TRUE)
tissue <- args[1] # do this for all tissues in GTEx

outdir='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex-v8/'
cov_base='/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/'

exprdir='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/'
eqtldir='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/'

################################################################
## filter GTEx VCF by SNPs in snp_list (snp_list == all significant associations in all tissues)
################################################################

vcf='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz'
snp_path <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/gtex.v8.esnp.txt' # all significant associations in GTEx v8
out <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/gtexv8_filtered'
system('echo "Generating filtered VCF"')
cmd <- sprintf('vcftools --gzvcf %s --out %s --snps %s --recode', vcf, out, snp_path)
system(cmd)
system(sprintf('sed -i "/^##/d" %s','/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/gtexv8_filtered.recode.vcf'))

################################################################
## Get LA for each GTEx v8 eSNP
################################################################

outdir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/'
# for file in `ls /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL | grep "egenes"`; do zcat $file | cut -f12 | sed -e '1d' >> /mnt/lab_data/montgomery/nicolerg/gtex_esnp_list.txt; done
snp_path <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex.v8.esnp.txt'

snp_list <- fread(snp_path, header=FALSE, sep='\t')
snp_list <- snp_list[,V1]

gtex_subj <- fread('/users/nicolerg/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
gtex_subj <- gtex_subj[,.(SUBJID,RACE)]
gtex_subj <- gtex_subj[SUBJID %in% SUBJ_OF_INTEREST]
# RACE
# 1=Asian
# 2=Black
# 3=White
gtex_subj <- gtex_subj[as.numeric(RACE) %in% c(1,2,3)]
gtex_subj_ids <- gtex_subj[,SUBJID]

# filter local ancestry 
local_all <- list()
i <- 1
for(chrom in 1:22){
	print(sprintf('working on chr%s',chrom))
	local <- paste0(cov_base,'/chr',chrom,'.hg19.maf0.localcov.tsv.gz')
	local <- fread(cmd=sprintf('zcat %s',local), sep='\t', header=TRUE)
	local <- local[,.(SUBJID, SNP_ID, AFR, ASN, EUR)]
	subjs <- unique(local[,SUBJID])
	# filter SNPs
	local <- local[SNP_ID %in% snp_list]
	
	# add local covariates for missing subjects (EUR, AFR, ASN only)
	missing <- gtex_subj_ids[!gtex_subj_ids %in% subjs]
	snps <- unique(local[,SNP_ID])

	missing_cov <- list()
	j <- 1
	for (id in missing){
		#SUBJID  SNP_ID  AFR     ASN     EUR
		SUBJID <- id 
		SNP_ID <- snps

		race <- gtex_subj[SUBJID == id, RACE]
		if (race == 1){
			ASN=2
			AFR=0
			EUR=0
		} else if (race == 2){
			ASN=0
			AFR=2
			EUR=0
		} else if (race == 3){
			ASN=0
			EUR=2
			AFR=0
		}

		new = data.table(data.frame(SUBJID=SUBJID, SNP_ID=SNP_ID, AFR=AFR, ASN=ASN, EUR=EUR))
		missing_cov[[j]] <- new

		j <- j+1
	}

	missing_cov <- data.table(rbindlist(missing_cov))
	local <- data.table(rbind(local, missing_cov))

	local_all[[i]] <- local
	i <- i+1

}

local_all_gtex <- data.table(rbindlist(local_all))
save(local_all_gtex, file=sprintf('%s/local_all_gtex_lead_snps.RData',outdir))

################################################################

# load local ancestry calls 
writeLines('Reading in local ancestry')
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/local_all_gtex_lead_snps.RData')

writeLines('Reading in VCF')
# read in VCF 
vcf <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/gtexv8_filtered.recode.vcf',sep='\t',header=T)
colnames(vcf)[3] <- 'SNP_ID'

# > head(local_all_gtex)
#        SUBJID             SNP_ID AFR ASN EUR
# 1: GTEX-11DYG chr1_58814_G_A_b38   1   0   0
# 2: GTEX-11NV4 chr1_58814_G_A_b38   0   0   0

# tissues <- system(sprintf('ls %s | sed "s/\\.v8.*//" | uniq', exprdir),intern=T)

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

writeLines('Filtering subjects')
# include only people in this tissue 
header <- fread(cmd=sprintf('zcat %s/%s.v8.normalized_expression.bed.gz',exprdir,tissue), sep='\t', header=T, nrows=1)
subjs <- colnames(header)[5:ncol(header)]
sub_local <- local_all_gtex[SUBJID %in% subjs]
local_all_gtex <- NULL # clear some memory

writeLines('Filtering SNPs')
# only include SNPs in this tissue
sig_assoc <- fread(cmd=sprintf('zcat %s/%s.v8.signif_variant_gene_pairs.txt.gz', eqtldir, tissue), sep='\t', header=T)
test_snps <- unique(sig_assoc[,variant_id])
sub_local <- sub_local[SNP_ID %in% test_snps]
sig_assoc <- NULL # clear some memory
print(length(test_snps))

writeLines('Reformatting VCF')
# subset VCF (both people and SNPs); reformat
vcf_filt <- vcf[SNP_ID %in% test_snps]
vcf <- NULL # clear some memory
cols <- subjs[subjs %in% colnames(vcf_filt)]
vcf_filt <- vcf_filt[,c('SNP_ID', cols),with=FALSE]
writeLines('Converting haplotypes to number')
vcf_filt[ , (cols) := lapply(.SD, convert_col), .SDcols = cols]
row_label <- vcf_filt[,SNP_ID]
vcf_filt[,SNP_ID := NULL]
subjid <- colnames(vcf_filt)
vcf_filt <- data.table(t(vcf_filt))
colnames(vcf_filt) <- row_label
vcf_filt[,SUBJID := subjid]

writeLines('Running regressions')

# create a progress bar
total <- length(test_snps)
count = 0
pb <- txtProgressBar(min = 0, max = total, style = 3)

outfile <- file(paste0(outdir,'/',tissue,'_leadsnp.r2.txt'), 'w')
cat("variant_id\tlocal_rsq\n", file=outfile)

try_lm <- function(f, df) {
    out <- tryCatch(
        {
			lm(f, data=df)$r.squared
        },
        error=function(cond) {
            return(NA)
        }
    )    
    return(out)
}

for (v in test_snps){
	
	# skip SNPs not in this file
	if (!v %in% colnames(vcf_filt)){
		next
	}

	g <- vcf_filt[,c('SUBJID', v), with=FALSE]
	l <- sub_local[SNP_ID == v]
	m <- merge(g, l, by='SUBJID')
	f <- as.formula(sprintf('%s ~ AFR + ASN',v))

	rsq <- try_lm(f, m)

	string <- paste(v, rsq, sep='\t')
	print(string)
	cat(paste0(string,'\n'), file=outfile, append=TRUE)

	count <- count + 1
	setTxtProgressBar(pb, count)
}
close(pb)
close(outfile)
