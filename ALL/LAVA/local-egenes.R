# find gene-variant associations for all pairs in the eGenes file (log10(p_val) > 5)

library(data.table)

genofile="/mnt/lab_data/montgomery/nicolerg/filtered/gtex.phased.MAF01.hg19.chr21.vcf"
globalcovfile="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
localcovfile="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/chr21.hg19.localcov.tsv"
phenofile="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
exprfile="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/Whole_Blood_expression_chr21.bed"
egenesfile="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/egenes.chr21.tsv"
out1="/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood-local-chr21-egenes.tsv"

system("echo 'Loading files'")
geno <- fread(genofile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
expr <- fread(exprfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
globalcov <- read.delim(globalcovfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
localcov <- fread(localcovfile, header=TRUE, sep='\t', stringsAsFactors=FALSE)
localcov <- data.frame(localcov, stringsAsFactors=FALSE)
pheno <- read.delim(file=phenofile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
egenes <- read.delim(file=egenesfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)

colnames(globalcov) <- gsub('.','-',colnames(globalcov),fixed=TRUE)
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

system("echo 'Sorting subjects in all tables'")
# select only White, African American, and Asian individuals from the VCF
pheno <- pheno[grepl("1|2|3",pheno$RACE),]
# also keep track of White
white <- pheno$SUBJID[grepl("3",pheno$RACE)]
m <- match(pheno$SUBJID, colnames(expr))
m <- na.omit(m)
names <- colnames(expr)[m]
expr <- expr[,c(1:4,m),with=FALSE]
# prune by those present in the expression file
m <- match(names, colnames(geno))
geno <- geno[,c(1:9,m),with=FALSE]
m <- match(names, colnames(globalcov))
globalcov <- globalcov[,c(1,m)]
white <- white[na.omit(match(names, white))]
localcov$EUR <- NULL

# now all df have only black, white, asian subjects also in expr, in the same order everywhere
geno <- data.frame(geno, stringsAsFactors=FALSE) # convert from data table to data frame
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

expr <- data.frame(expr, stringsAsFactors=FALSE)
colnames(expr) <- gsub('.','-',colnames(expr),fixed=TRUE)

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

maf <- rowMeans(geno[,10:ncol(geno)])/2
maf <- pmin(maf, 1-maf)

# transpose global covariates
rownames(globalcov) <- globalcov$ID
globalcov$ID <- NULL
globalcov <- t(globalcov)

globalcov <- globalcov[,-c(1:5)]

# make a df for "local covariates" for Whites
system("echo")
localwhite <- data.frame(matrix(ncol=3,nrow=length(white)))
rownames(localwhite) <- white
colnames(localwhite) <- colnames(localcov)[4:ncol(localcov)]
localwhite$AFR <- 0
localwhite$ASN <- 0
localwhite$UNK <- 0
allpairs <- data.frame(matrix(ncol=15, nrow=0))
colnames(allpairs) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","afr_pval_nominal","afr_slope","afr_slope_se","asn_pval_nominal","asn_slope","asn_slope_se")
# initialize data frame to keep track of "unknowns" removed for each 
unknown_count <- data.frame(matrix(ncol=2, nrow=0))
colnames(unknown_count) <- c("SNP_ID","FRAC_REMOVED")
system("echo 'Start eQTL calling with global and local covariates'")
# remove PCs from global covariates
globalcov <- globalcov[,-c(1:5)]

for(j in 1:nrow(egenes)){
	gene_id <- egenes$gene_id[j]
	variant_id <- egenes$variant_id[j]
	print(variant_id)
	local <- localcov[grepl(variant_id,localcov$SNP_ID),] # just admixed
	n <- local$SUBJID
	local <- local[,4:ncol(local)]
	rownames(local) <- n
	alocal <- rbind(local, localwhite)
	n <- c(rownames(local),rownames(localwhite))
	rownames(alocal) <- n
	# select individuals only in globalcov
	m <- match(rownames(globalcov),rownames(alocal))
	n <- rownames(alocal)[m]

	if((length(n) - length(na.omit(n))) > 0){

		print(paste("Variant",variant_id, "did not convert to hg19",sep=' '))
		row <- data.frame(matrix(c(gene_id, variant_id, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'),nrow=1,ncol=15))
		colnames(row) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","afr_pval_nominal","afr_slope","afr_slope_se","asn_pval_nominal","asn_slope","asn_slope_se")
		allpairs <- rbind(allpairs, row)

	} else {

		alocal <- alocal[m,]
		alocal <- data.frame(lapply(alocal, as.numeric),stringsAsFactors=FALSE)
		rownames(alocal) <- n

		# extract genotype data for variant_id
		row <- grep(variant_id, geno$ID, fixed=TRUE)
		if(length(row)==0){
			print(paste("No genotype data for",variant_id,sep=" "))
			row <- data.frame(matrix(c(gene_id, variant_id, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'),nrow=1,ncol=15))
			colnames(row) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","afr_pval_nominal","afr_slope","afr_slope_se","asn_pval_nominal","asn_slope","asn_slope_se")
			allpairs <- rbind(allpairs, row)

		} else {

			m <- unname(maf[row])
			ma_count = sum(geno[row,][-c(1:9)])
			ma_samples = length(grep("1|2",geno[row,10:ncol(geno)]))
			haplo <- data.frame(unlist(geno[row,])[-c(1:9)],stringsAsFactors=FALSE)
			colnames(haplo) <- "haplo"
			varpos <- as.numeric(unlist(strsplit(geno$ID[row], '_'))[2])

			# extract expression data for gene_id
			row <- grep(gene_id, expr$gene_id, fixed=TRUE)

			if(length(row)==0){

				print(paste("No expression data for",gene_id,sep=" "))
				row <- data.frame(matrix(c(gene_id, variant_id, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'),nrow=1,ncol=15))
				colnames(row) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","afr_pval_nominal","afr_slope","afr_slope_se","asn_pval_nominal","asn_slope","asn_slope_se")
				allpairs <- rbind(allpairs, row)

			} else {

				expression = data.frame(as.numeric(unlist(expr[row,])[-c(1:4)]),stringsAsFactors=FALSE)
				colnames(expression) <- "expression"
				tss_distance = as.numeric(varpos - as.numeric(expr$start[row]))

				# put all the parts togther
				r <- rownames(globalcov)
				df1 <- data.frame(cbind(expression, haplo, globalcov),stringsAsFactors=FALSE)
				rownames(df1) <- r
				df1 <- df1[match(rownames(alocal),rownames(df1)),]
				r <- rownames(alocal)
				df <- data.frame(cbind(df1, alocal),stringsAsFactors=FALSE)
				df[,1:62] <- data.frame(lapply(df[,1:62],as.numeric),stringsAsFactors=FALSE)

				df$ASN <- as.numeric(df$ASN)
				df$AFR <- as.numeric(df$AFR)

				df$platform <- as.factor(df$platform)
				df$pcr <- as.factor(df$pcr)
				df$sex <- as.factor(df$sex)
				# # remove factor columns that have only one level
				# remove <- c()
				# for(z in (ncol(df)-5):(ncol(df))){
				# 	if(length(levels(df[,z])) == 1){
				# 		remove <- c(remove, z)
				# 	}
				# }
				# if(length(remove)>0){
				# 	df <- df[,-remove]
				# }
				rownames(df) <- r
				# linear regression
				lm.fit <- lm(expression ~ ., data=df)
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

				print(afr_pval_nominal)
				
				row <- data.frame(matrix(c(gene_id, variant_id, tss_distance, ma_samples, ma_count, m, pval_nominal, slope, slope_se, afr_pval_nominal, afr_slope, afr_slope_se, asn_pval_nominal, asn_slope, asn_slope_se),nrow=1,ncol=15))
				colnames(row) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","afr_pval_nominal","afr_slope","afr_slope_se","asn_pval_nominal","asn_slope","asn_slope_se")
				allpairs <- rbind(allpairs, row)
			}

		}

	}
}

system("echo 'Writing all_pairs file for eQTL calling with global and local covariates'")
write.table(allpairs, out1, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
