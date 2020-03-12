library(data.table)

args <- commandArgs(trailingOnly=TRUE)
genofile <- args[1]
exprfile <- args[2]
globalcovfile <- args[3]
localcovfile <- args[4]
phenofile <- args[5]
outfile <- args[6]

# genofile <- "/mnt/lab_data/montgomery/nicolerg/filtered/gtex.phased.MAF01.hg19.chr21.vcf"
# exprfile <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/Whole_Blood_expression_chr21.bed"
# globalcovfile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
# localcovfile <- "/mnt/lab_data/montgomery/nicolerg/admix/chr21_subsets/chr21.hg19.localcov.tsv"
# phenofile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

system("echo 'Loading files'")
geno <- fread(genofile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
expr <- fread(exprfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
globalcov <- read.delim(globalcovfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
localcov <- fread(localcovfile, header=TRUE, sep='\t', stringsAsFactors=FALSE)
pheno <- read.delim(file=phenofile, sep='\t', header=TRUE, stringsAsFactors=FALSE)

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
localnames <- unique(localcov$SUBJID)
for(i in 1:length(localnames)){
	name = localnames[i]
	if(length(grep(name, names))==0){
		localcov <- localcov[!grepl(name, localcov$SUBJID),]
	}
}
localcov$EUR <- NULL

# now all df have only black, white, asian subjects also in expr, in the same order everywhere
geno <- data.frame(geno, stringsAsFactors=FALSE) # convert from data table to data frame
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

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

# remove variants that are too far from this subset of genes
system("echo 'Removing variants that are too far from the subset of genes'")
minstart <- min(as.numeric(expr$start)) - 1500000 # 500 kb buffer to accommodate different genome builds
maxstart <- max(as.numeric(expr$start)) + 1500000
geno <- geno[!((geno$POS < minstart) | (geno$POS > maxstart)),]
expr <- data.frame(expr, stringsAsFactors=FALSE)
colnames(expr) <- gsub('.','-',colnames(expr),fixed=TRUE)
maf <- rowMeans(geno[,10:ncol(geno)])/2
maf <- pmin(maf, 1-maf)
# # remove genes that are too far from this sample of variants
# system("echo 'Removing genes that are too far from the subset of SNPs'")
# minstart <- min(as.numeric(geno$POS)) - 1500000
# maxstart <- max(as.numeric(geno$POS)) + 1500000
# expr <- expr[!((expr$start < minstart) | (expr$start > maxstart)),] # 24 genes 
# expr <- data.frame(expr, stringsAsFactors=FALSE)
# colnames(expr) <- gsub('.','-',colnames(expr),fixed=TRUE)
# maf <- rowMeans(geno[,10:ncol(geno)])/2
# maf <- pmin(maf, 1-maf)
# # geno <- geno[(maf > 0.01),] # 5222 variants
# # print(nrow(geno))
# # maf <- rowMeans(geno[,10:ncol(geno)])/2
# # maf <- pmin(maf, 1-maf)
# # transpose global covariates
rownames(globalcov) <- globalcov$ID
globalcov$ID <- NULL
globalcov <- t(globalcov)

globalcov <- globalcov[,-c(1:5)]
# gid <- data.frame(rownames(globalcov),pheno$RACE[match(rownames(globalcov),pheno$SUBJID)])
# write.table(gid,"global-ids.txt",sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # eQTL calling with no local covariates
# allpairs <- data.frame(matrix(ncol=9, nrow=0))
# colnames(allpairs) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se")
# system("echo 'Start eQTL calling with only global covariates'")
# for(i in 1:nrow(expr)){
# 	for(j in 1:nrow(geno)){
# 		# check that the SNP is within 1 MB of the gene TSS 
# 		varpos <- as.numeric(unlist(strsplit(geno$ID[j], '_'))[2])
# 		tss_distance = as.numeric(varpos - as.numeric(expr$start[i]))
# 		if(tss_distance <= 1000000 & tss_distance >= -1000000){
# 			# limit to ma_samples >= 10 
# 			ma_samples = length(grep("1|2",geno[j,10:ncol(geno)]))
# 			if(ma_samples >= 10){
# 				# define statistics for output file
# 				variant_id = geno$ID[j]
# 				gene_id = expr$gene_id[i]
# 				print(paste("Working on gene",gene_id,"and SNP",variant_id,"..."),sep=' ')
# 				ma_count = sum(geno[j,][-c(1:9)])
# 				# generate data frame for regression
# 				expression = data.frame(unlist(expr[i,])[-c(1:4)],stringsAsFactors=FALSE)
# 				colnames(expression) <- "expression"
# 				haplo = data.frame(unlist(geno[j,10:ncol(geno)]),stringsAsFactors=FALSE)
# 				colnames(haplo) <- "haplo"
# 				r <- rownames(globalcov)
# 				df <- data.frame(cbind(expression, haplo, globalcov),stringsAsFactors=FALSE)
# 				df <- data.frame(lapply(df, as.numeric))
# 				df$platform <- as.factor(df$platform)
# 				df$pcr <- as.factor(df$pcr)
# 				df$sex <- as.factor(df$sex)
# 				rownames(df) <- r
# 				# linear regression
# 				lm.fit <- lm(expression ~ ., data=df)
# 				if(!is.null(warnings())){
# 					print(warnings())
# 				}
# 				out <- summary(lm.fit)$coefficients
# 				pval_nominal = out["haplo",4]
# 				slope = out["haplo",1]
# 				slope_se = out["haplo",2]
# 				m <- unname(maf[j])
# 				row <- data.frame(matrix(c(gene_id, variant_id, tss_distance, ma_samples, ma_count, m, pval_nominal, slope, slope_se),nrow=1,ncol=9))
# 				colnames(row) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se")
# 				allpairs <- rbind(allpairs, row)
# 			}
# 		}
# 	}
# }
# system("echo 'Writing all_pairs file for eQTL calling with global covariates'")
# save(allpairs, file="/mnt/lab_data/montgomery/nicolerg/admix/wholeblood-chr21subset-allpairs-globalcov-noPCs.RData")
# write.table(allpairs, "/mnt/lab_data/montgomery/nicolerg/admix/wholeblood-chr21subset-allpairs-globalcov-hg19-500-noPCs.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

# colnames(globalcov)[1:20] # col 1: ID
# colnames(expr)[1:20] # col 1-4: #chr, start, end, gene_id
# colnames(geno)[1:20] # col 1-9: POS, ID

# # now, do it with local covariates
# # make a df for "local covariates" for Whites
# system("echo")
# #localwhite <- data.frame(matrix(ncol=8,nrow=length(white)))
# localwhite <- data.frame(matrix(ncol=3,nrow=length(white)))
# rownames(localwhite) <- white
# colnames(localwhite) <- colnames(localcov)[4:ncol(localcov)]
# localwhite$AFR <- 0
# localwhite$ASN <- 0
# localwhite$UNK <- 0
# # localwhite$EUR <- 2
# allpairs <- data.frame(matrix(ncol=15, nrow=0))
# colnames(allpairs) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","afr_pval_nominal","afr_slope","afr_slope_se","asn_pval_nominal","asn_slope","asn_slope_se")
# # initialize data frame to keep track of "unknowns" removed for each 
# unknown_count <- data.frame(matrix(ncol=2, nrow=0))
# colnames(unknown_count) <- c("SNP_ID","FRAC_REMOVED")
# system("echo 'Start eQTL calling with global and local covariates'")
# # remove PCs from global covariates
# globalcov <- globalcov[,-c(1:5)]
# for(j in 1:nrow(geno)){
# 	# make data table for local covariates
# 	variant_id = geno$ID[j]
# 	local <- localcov[grepl(variant_id,localcov$SNP_ID),] # just admixed
# 	n <- local$SUBJID
# 	local <- local[,4:ncol(local)]
# 	rownames(local) <- n
# 	alocal <- rbind(local, localwhite)
# 	n <- c(rownames(local),rownames(localwhite))
# 	rownames(alocal) <- n
# 	# select individuals only in globalcov
# 	m <- match(rownames(globalcov),rownames(alocal))
# 	alocal <- alocal[m,]
# 	alocal <- data.frame(lapply(alocal, as.numeric),stringsAsFactors=FALSE)
# 	rownames(alocal) <- n
# 	print(nrow(alocal))
# 	# remove rows with UNK; keep track of number
# 	ninit <- nrow(alocal)
# 	alocal <- alocal[grep(0,alocal$UNK),]
# 	nfin <- nrow(alocal)
# 	diff <- (ninit-nfin)/ninit
# 	alocal$UNK <- NULL
# 	row <- data.frame(matrix(c(variant_id,diff),nrow=1,ncol=2))
# 	colnames(row) <- c("SNP_ID","FRAC_REMOVED")
# 	unknown_count <- rbind(unknown_count, row)
# 	for(i in 1:nrow(expr)){
# 		# check that the SNP is within 1 MB of the gene TSS 
# 		varpos <- as.numeric(unlist(strsplit(geno$ID[j], '_'))[2])
# 		tss_distance = as.numeric(varpos - as.numeric(expr$start[i]))
# 		if(tss_distance <= 1000000 & tss_distance >= -1000000){
# 			# limit to ma_samples >= 10 
# 			ma_samples = length(grep("1|2",geno[j,10:ncol(geno)]))
# 			if(ma_samples >= 10){
# 				# define statistics for output file
# 				variant_id = geno$ID[j]
# 				gene_id = expr$gene_id[i]
# 				print(paste("Working on gene",gene_id,"and SNP",variant_id,"..."),sep=' ')
# 				ma_count = sum(geno[j,][-c(1:9)])
# 				# generate data frame for regression
# 				expression = data.frame(as.numeric(unlist(expr[i,])[-c(1:4)]),stringsAsFactors=FALSE)
# 				colnames(expression) <- "expression"
# 				haplo = data.frame(unlist(geno[j,])[-c(1:9)],stringsAsFactors=FALSE)
# 				colnames(haplo) <- "haplo"
# 				# put all the parts togther
# 				r <- rownames(globalcov)
# 				df1 <- data.frame(cbind(expression, haplo, globalcov),stringsAsFactors=FALSE)
# 				rownames(df1) <- r
# 				df1 <- df1[match(rownames(alocal),rownames(df1)),]
# 				r <- rownames(alocal)
# 				df <- data.frame(cbind(df1, alocal),stringsAsFactors=FALSE)
# 				df[,1:62] <- data.frame(lapply(df[,1:62],as.numeric),stringsAsFactors=FALSE)

# 				df$ASN <- as.numeric(df$ASN)
# 				df$AFR <- as.numeric(df$AFR)

# 				df$platform <- as.factor(df$platform)
# 				df$pcr <- as.factor(df$pcr)
# 				df$sex <- as.factor(df$sex)
# 				# remove factor columns that have only one level
# 				remove <- c()
# 				#for(z in (ncol(df)-6):(ncol(df))){
# 				for(z in (ncol(df)-10):(ncol(df))){
# 					if(length(levels(df[,z])) == 1){
# 						remove <- c(remove, z)
# 					}
# 				}
# 				if(length(remove)>0){
# 					df <- df[,-remove]
# 				}
# 				rownames(df) <- r
# 				# linear regression
# 				lm.fit <- lm(expression ~ ., data=df)
# 				#lm.fit <- lm(expression ~ . + haplo:ASN + haplo:AFR, data=df)
# 				if(!is.null(warnings())){
# 					print(warnings())
# 				}
# 				out <- summary(lm.fit)$coefficients
				
# 				pval_nominal = out["haplo",4]
# 				slope = out["haplo",1]
# 				slope_se = out["haplo",2]

# 				afr_pval_nominal = out["AFR",4]
# 				afr_slope = out["AFR",1]
# 				afr_slope_se = out["AFR",2]

# 				asn_pval_nominal = out["ASN",4]
# 				asn_slope = out["ASN",1]
# 				asn_slope_se = out["ASN",2]

# 				m <- unname(maf[j])
# 				row <- data.frame(matrix(c(gene_id, variant_id, tss_distance, ma_samples, ma_count, m, pval_nominal, slope, slope_se, afr_pval_nominal, afr_slope, afr_slope_se, asn_pval_nominal, asn_slope, asn_slope_se),nrow=1,ncol=15))
# 				colnames(row) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","afr_pval_nominal","afr_slope","afr_slope_se","asn_pval_nominal","asn_slope","asn_slope_se")
# 				allpairs <- rbind(allpairs, row)
# 			}
# 		}
# 	}
# }
# save(allpairs, file="/mnt/lab_data/montgomery/nicolerg/admix/wholeblood-chr21subset-allpairs-localcov.RData")
# system("echo 'Writing all_pairs file for eQTL calling with global and local covariates'")
# write.table(allpairs, "/mnt/lab_data/montgomery/nicolerg/admix/wholeblood-chr21subset-allpairs-localcov-hg19-500-noPCs.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
# write.table(unknown_count, "/mnt/lab_data/montgomery/nicolerg/admix/wholeblood-chr21subset-unknown_count.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
