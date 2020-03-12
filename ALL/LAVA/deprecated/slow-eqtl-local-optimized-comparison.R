library(data.table)

args <- commandArgs(trailingOnly=TRUE)
genofile <- args[1]
exprfile <- args[2]
globalcovfile <- args[3]
localcovfile <- args[4]
phenofile <- args[5]
out1 <- args[6]

# genofile <- "/mnt/lab_data/montgomery/nicolerg/filtered/gtex.phased.MAF01.hg19.chr21.vcf"
# exprfile <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/Whole_Blood_expression_chr21.bed"
# globalcovfile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
# localcovfile <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/chr21.hg19.localcov.tsv"
# phenofile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

system("echo 'Loading files'")
geno <- fread(input=sprintf("zcat %s", genofile), sep='\t', header=TRUE, stringsAsFactors=FALSE)
expr <- fread(input=sprintf("zcat %s", exprfile), sep='\t', header=TRUE, stringsAsFactors=FALSE)
globalcov <- read.delim(globalcovfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
localcov <- fread(input=sprintf("zcat %s", localcovfile), header=TRUE, sep='\t', stringsAsFactors=FALSE)
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
colnames(globalcov)

# make a df for "local covariates" for Whites
system("echo")
localwhite <- data.table(SUBJID=white, AFR=rep(0, length(white)), ASN=rep(0, length(white)), UNK=rep(0, length(white)))

globalcov <- data.table(globalcov, keep.rownames="SUBJID")
setkey(globalcov, SUBJID)

allp <- file(out1, 'w')
cat("gene_id\tvariant_id\tglobal_pval\tglobal_slope\tglobal_slope_se\tpc1_pval\tpc2_pval\tpc3_pval\tpc4_pval\tpc5_pval\tlocal_pval\tlocal_slope\tlocal_slope_se\tafr_pval\tafr_slope\tafr_slope_se\tasn_pval\tasn_slope\tasn_slope_se\tanova_pval\tanova_global_rss\tanova_local_rss\n", file=allp)

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

				df <- data.frame(dt, stringsAsFactors=FALSE)

				# regression without local covariates
				global.fit <- lm(expression ~ . - AFR - ASN - SUBJID, data=df)

				if(!is.null(warnings())){
					print(warnings())
				}
				out <- summary(global.fit)$coefficients

				pc1_pval = out["PC1",4]
				pc2_pval = out["PC2",4]
				pc3_pval = out["PC3",4]
				pc4_pval = out["PC4",4]
				pc5_pval = out["PC5",4]

				if ("haplo"%in%rownames(out)==FALSE) {
					global_pval = "NA"
					gloval_slope = "NA"
					global_slope_se = "NA"
				} else {
					global_pval = out["haplo",4]
					global_slope = out["haplo",1]
					global_slope_se = out["haplo",2]
				}

				# regression with local covariates
				local.fit <- lm(expression ~ . - SUBJID, data=df)

				if(!is.null(warnings())){
					print(warnings())
				}
				out <- summary(local.fit)$coefficients

				if ("ASN"%in%rownames(out)==FALSE) {
					asn_pval = "NA"
					asn_slope = "NA"
					asn_slope_se = "NA"
				} else {
					asn_pval = out["ASN",4]
					asn_slope = out["ASN",1]
					asn_slope_se = out["ASN",2]
				}

				if ("AFR"%in%rownames(out)==FALSE) {
					afr_pval = "NA"
					afr_slope = "NA"
					afr_slope_se = "NA"
				} else {
					afr_pval = out["AFR",4]
					afr_slope = out["AFR",1]
					afr_slope_se = out["AFR",2]
				}

				if ("haplo"%in%rownames(out)==FALSE) {
					local_pval = "NA"
					local_slope = "NA"
					local_slope_se = "NA"
				} else {
					local_pval = out["haplo",4]
					local_slope = out["haplo",1]
					local_slope_se = out["haplo",2]
				}

				# compare models with anova
				a <- anova(global.fit,local.fit)
				anova_pval <- unname(unlist(a["Pr(>F)"][2,]))
				anova_local_rss <- unname(unlist(a["RSS"][2,]))
				anova_global_rss <- unname(unlist(a["RSS"][1,]))

				string=paste(gene_id, variant_id, global_pval, global_slope, global_slope_se, pc1_pval, pc2_pval, pc3_pval, pc4_pval, pc5_pval, local_pval, local_slope, local_slope_se, afr_pval, afr_slope, afr_slope_se, asn_pval, asn_slope, asn_slope_se, anova_pval, anova_global_rss, anova_local_rss, sep='\t')
				cat(paste(string,'\n',sep=''), file=allp, append=TRUE)
			}
		}
	} else {
		count = count + nrow(expr)
	}
}
close(pb)
close(allp)
