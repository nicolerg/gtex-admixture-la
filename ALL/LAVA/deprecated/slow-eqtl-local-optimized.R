library(data.table)

args <- commandArgs(trailingOnly=TRUE)
genofile <- args[1]
exprfile <- args[2]
globalcovfile <- args[3]
localcovfile <- args[4]
phenofile <- args[5]
out1 <- args[6]
out2 <- args[7]

# genofile <- "/mnt/lab_data/montgomery/nicolerg/filtered/gtex.phased.MAF01.hg19.chr21.vcf"
# exprfile <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/Whole_Blood_expression_chr21.bed"
# globalcovfile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
# localcovfile <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/chr21.hg19.localcov.tsv"
# phenofile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

system("echo 'Loading files'")
geno <- fread(input=sprintf("zcat %s", genofile), sep='\t', header=TRUE, stringsAsFactors=FALSE)
expr <- fread(exprfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
globalcov <- read.delim(globalcovfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
localcov <- fread(localcovfile, header=TRUE, sep='\t', stringsAsFactors=FALSE)
pheno <- fread(input=sprintf("cat %s", phenofile), sep='\t', header=TRUE)

colnames(globalcov) <- gsub('.','-',colnames(globalcov),fixed=TRUE)
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

system("echo 'Sorting subjects in all tables'")
# select only White, African American, and Asian individuals from the VCF (pheno data table is used only to filter SUBJID)
pheno <- pheno[ RACE == '1' | RACE == '2' | RACE == '3' ]
# also keep track of White
white <- pheno[ RACE == '3' , SUBJID ]
# filter expression file
m <- match(pheno[,SUBJID], colnames(expr))
m <- na.omit(m)
names <- colnames(expr)[m]
expr <- expr[,c(1:4,m),with=FALSE]
# filter genotype file
m <- match(names, colnames(geno))
geno <- geno[,c(1:9,m),with=FALSE]
m <- match(names, colnames(globalcov))
globalcov <- globalcov[,c(1,m)]
white <- white[na.omit(match(names, white))]
localcov <- localcov[, !c("EUR"), with = FALSE]

system("echo 'Filtering variants'")
# remove variants that are too far from this subset of genes
extract_id <- function(x){
	return(as.numeric(unlist(strsplit(x, '_'))[2]))
}
geno <- data.frame(geno, stringsAsFactors=FALSE, check.names=FALSE)
geno$POS <- sapply(geno$ID, extract_id)
minstart <- min(as.numeric(expr[,start])) - 1000000 
maxstart <- max(as.numeric(expr[,start])) + 1000000
geno <- data.table(geno)
# select by POS (data.table)
geno <- geno[as.integer(POS) >= minstart & as.integer(POS) <= maxstart]
geno <- data.frame(geno, stringsAsFactors=FALSE, check.names=FALSE)

system("echo 'Converting phased haplotypes to genotype dosage'")
# convert haplotypes to c(0,1,2)
hap2bin <- function(x){
	if (x == '0|0'){
		return(0)
	} else if (x == '1|1'){
		return(2)
	} else {
		return(1)
	}
}
vec_h2b <- function(xs){
	sapply(xs, hap2bin)
}
for (j in 10:ncol(geno)){
	new <- vec_h2b(geno[,j])
	geno[,j] <- new
}

# remove variants with minor allele in less than 10 samples
geno <- data.table(geno)
geno <- geno[rowSums(geno[,10:ncol(geno)]) >= 10]
maf <- rowMeans(geno[,10:ncol(geno),with=FALSE])/2
maf <- pmin(maf, 1-maf)
print(paste('N filtered variants:',nrow(geno),sep=' '))

# transpose global covariates
rownames(globalcov) <- globalcov$ID
globalcov$ID <- NULL
globalcov <- t(globalcov)

# make a df for "local covariates" for Whites
system("echo")
localwhite <- data.table(SUBJID=white, AFR=rep(0, length(white)), ASN=rep(0, length(white)), UNK=rep(0, length(white)))

# remove PCs from global covariates
globalcov <- globalcov[,-c(1:5)]
globalcov <- data.table(globalcov, keep.rownames="SUBJID")
setkey(globalcov, SUBJID)

allp <- file(out1, 'w')
unkcounts <- file(out2, 'w')
cat("gene_id\tvariant_id\ttss_distance\tma_samples\tma_count\tmaf\tpval_nominal\tslope\tslope_se\tafr_pval_nominal\tafr_slope\tafr_slope_se\tasn_pval_nominal\tasn_slope\tasn_slope_se\n", file=allp)
cat("SNP_ID\tFRAC_REMOVED\n", file=unkcounts)

# data tables: expr, geno, localwhite, localcov, globalcov
# data frames: 

# create a progress bar
total <- nrow(geno)*nrow(expr)
count = 0
system("echo 'Start eQTL calling with global and local covariates'")
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(j in 1:nrow(geno)){
	setTxtProgressBar(pb, count)
	ma_samples = length(grep('1|2', geno[j]))
	if(ma_samples >= 10){
		# make data table for local covariates
		variant_id = geno[j,ID]
		local <- localcov[SNP_ID == variant_id] # this should be the biggest time saver
		local <- local[,.(SUBJID,AFR,ASN,UNK)]
		l = list(local, localwhite)
		alocal <- rbindlist(l , use.names=TRUE)
		setkey(alocal, SUBJID)
		strTmp = c('AFR','ASN','UNK') # Name of columns to be converted to numeric
		alocal[, (strTmp) := lapply(.SD, as.integer), .SDcols = strTmp]
		# remove rows with UNK; keep track of number
		ninit <- nrow(alocal)
		all_local <- alocal[UNK == 0]
		nfin <- nrow(all_local)
		diff <- (ninit-nfin)/ninit
		all_local <- all_local[, !c('UNK'), with = FALSE]
		setkey(all_local, SUBJID)
		# keep track of how many are removed
		s = paste(variant_id, diff, sep='\t')
		cat(paste(s,'\n', sep=''), file=unkcounts, append=TRUE)
		for(i in 1:nrow(expr)){
			count = count + 1
			# check that the SNP is within 1 MB of the gene TSS 
			pos <- as.numeric(geno[j,POS])
			tss_distance = as.numeric(pos - expr[i,start])
			if(tss_distance <= 1000000 & tss_distance >= -1000000){
				# define statistics for output file
				gene_id = expr[i,gene_id]
				# generate data frame for regression
				df1 <- data.table(expression=as.numeric(unname(unlist(expr[i]))[-c(1:4)]),haplo=as.numeric(unname(unlist(geno[j]))[-c(1:9)]))
				df2 <- cbind(df1, globalcov)
				setkey(df2, SUBJID)
				dt <- df2[all_local,nomatch=0]
				strTmp = c('platform','sex','pcr')
				dt[, (strTmp) := lapply(.SD, as.factor), .SDcols = strTmp]
				ma_counts = sum(dt[,haplo])
				# # remove factor columns that have only one level
				# remove <- c()
				# #for(z in (ncol(df)-6):(ncol(df))){
				# for(z in (ncol(df)-10):(ncol(df))){
				# 	if(length(levels(df[,z])) == 1){
				# 		remove <- c(remove, z)
				# 	}
				# }
				# if(length(remove)>0){
				# 	df <- df[,-remove]
				# }
				# linear regression
				df <- data.frame(dt, stringsAsFactors=FALSE)

				lm.fit <- lm(expression ~ . - SUBJID, data=df)

				#lm.fit <- lm(expression ~ . + haplo:ASN + haplo:AFR, data=df)
				if(!is.null(warnings())){
					print(warnings())
				}
				out <- summary(lm.fit)$coefficients

				if ("ASN"%in%rownames(out)==FALSE) {
					asn_pval_nominal = "NA"
					asn_slope = "NA"
					asn_slope_se = "NA"
				} else {
					asn_pval_nominal = out["ASN",4]
					asn_slope = out["ASN",1]
					asn_slope_se = out["ASN",2]
				}

				if ("AFR"%in%rownames(out)==FALSE) {
					afr_pval_nominal = "NA"
					afr_slope = "NA"
					afr_slope_se = "NA"
				} else {
					afr_pval_nominal = out["AFR",4]
					afr_slope = out["AFR",1]
					afr_slope_se = out["AFR",2]
				}

				if ("haplo"%in%rownames(out)==FALSE) {
					pval_nominal = "NA"
					slope = "NA"
					slope_se = "NA"
				} else {
					pval_nominal = out["haplo",4]
					slope = out["haplo",1]
					slope_se = out["haplo",2]
				}

				m <- unname(maf[j])
				string=paste(gene_id, variant_id, tss_distance, ma_samples, ma_counts, m, pval_nominal, slope, slope_se, afr_pval_nominal, afr_slope, afr_slope_se, asn_pval_nominal, asn_slope, asn_slope_se, sep='\t')
				cat(paste(string,'\n',sep=''), file=allp, append=TRUE)
			}
		}
	} else {
		count = count + nrow(expr)
	}
}
close(pb)
close(allp)
close(unkcounts)
# system("echo 'Writing all_pairs file for eQTL calling with global and local covariates'")
# write.table(allpairs, out1, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
# write.table(unknown_count, out2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
