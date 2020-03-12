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

# now all df have only black, white, asian subjects also in expr, in the same order everywhere
geno <- data.frame(geno, stringsAsFactors=FALSE) # convert from data table to data frame
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

# remove variants that are too far from this subset of genes
system("echo 'Removing variants that are too far from the subset of genes'")
minstart <- min(as.numeric(expr$start)) - 1500000 # 500 kb buffer to accommodate different genome builds
maxstart <- max(as.numeric(expr$start)) + 1500000
print(nrow(geno))
geno <- geno[!((geno$POS < minstart) | (geno$POS > maxstart)),]
print(nrow(geno))
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

# eQTL calling with no local covariates
allpairs <- data.frame(matrix(ncol=9, nrow=0))
colnames(allpairs) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se")
system("echo 'Start eQTL calling with only global covariates'")
for(i in 1:nrow(expr)){
	for(j in 1:nrow(geno)){
		# check that the SNP is within 1 MB of the gene TSS 
		varpos <- as.numeric(unlist(strsplit(geno$ID[j], '_'))[2])
		tss_distance = as.numeric(varpos - as.numeric(expr$start[i]))
		if(tss_distance <= 1000000 & tss_distance >= -1000000){
			# limit to ma_samples >= 10 
			ma_samples = length(grep("1|2",geno[j,10:ncol(geno)]))
			if(ma_samples >= 10){
				# define statistics for output file
				variant_id = geno$ID[j]
				gene_id = expr$gene_id[i]
				#print(paste("Working on gene",gene_id,"and SNP",variant_id,"..."),sep=' ')
				ma_count = sum(geno[j,][-c(1:9)])
				# generate data frame for regression
				expression = data.frame(unlist(expr[i,])[-c(1:4)],stringsAsFactors=FALSE)
				colnames(expression) <- "expression"
				haplo = data.frame(unlist(geno[j,10:ncol(geno)]),stringsAsFactors=FALSE)
				colnames(haplo) <- "haplo"
				r <- rownames(globalcov)
				df <- data.frame(cbind(expression, haplo, globalcov),stringsAsFactors=FALSE)
				df <- data.frame(lapply(df, as.numeric))
				df$platform <- as.factor(df$platform)
				df$pcr <- as.factor(df$pcr)
				df$sex <- as.factor(df$sex)
				rownames(df) <- r
				# linear regression
				lm.fit <- lm(expression ~ ., data=df)
				if(!is.null(warnings())){
					print(warnings())
				}
				out <- summary(lm.fit)$coefficients
				if ("haplo"%in%rownames(out)==FALSE) {
					pval_nominal = "NA"
					slope = "NA"
					slope_se = "NA"
				} else {
					pval_nominal = out["haplo",4]
					slope = out["haplo",1]
					slope_se = out["haplo",2]
				}
				print(pval_nominal)
				m <- unname(maf[j])
				row <- data.frame(matrix(c(gene_id, variant_id, tss_distance, ma_samples, ma_count, m, pval_nominal, slope, slope_se),nrow=1,ncol=9))
				colnames(row) <- c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se")
				allpairs <- rbind(allpairs, row)
			}
		}
	}
}
system("echo 'Writing all_pairs file for eQTL calling with global covariates'")
write.table(allpairs, outfile, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
