library(data.table)

args <- commandArgs(trailingOnly=TRUE)
genofile <- args[1]
exprfile <- args[2]
globalcovfile <- args[3]
localcovfile <- args[4]
phenofile <- args[5]
out <- args[6]

system("echo 'Loading files'")
geno <- fread(input=sprintf("zcat %s", genofile), sep='\t', header=TRUE)
expr <- fread(input=sprintf("zcat %s", exprfile), sep='\t', header=TRUE)
globalcov <- read.delim(globalcovfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
localcov <- fread(input=sprintf("zcat %s", localcovfile), header=TRUE, sep='\t')
pheno <- fread(input=sprintf("cat %s", phenofile), sep='\t', header=TRUE)

colnames(globalcov) <- gsub('.','-',colnames(globalcov),fixed=TRUE)
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

system("echo 'Sorting subjects in all tables'")
# filter pheno by expr
subjs <- colnames(expr)[5:ncol(expr)]
pheno <- pheno[SUBJID %in% subjs]
# only keep AA individuals
# 1 = Asian
# 2 = Black
# 3 = White
nadx <- length(grep("2", pheno[,RACE]))
print(paste0("N black: ", length(grep("2", pheno[,RACE]))), quote=FALSE)
# filter expression file
subjs <- c(pheno[ RACE == '2', SUBJID ])
cols <- c("#chr","start","end","gene_id",subjs)
expr <- expr[,cols,with=FALSE]
# filter geno file
geno <- geno[,c(colnames(geno)[1:9], subjs), with=FALSE]
# filter global cov file
globalcov <- globalcov[,c("ID",subjs),with=FALSE]

system("echo 'Filtering variants'")
# remove variants that are too far from entire subset
extract_id <- function(x){
	return(as.numeric(unlist(strsplit(x, '_'))[2]))
}
minstart <- min(as.numeric(expr[,start])) - 1000000 
maxstart <- max(as.numeric(expr[,start])) + 1000000
geno <- geno[as.integer(POS) >= minstart & as.integer(POS) <= maxstart]

system("echo 'Converting phased haplotypes to genotype dosage'")
# convert haplotypes to c(0,1,2)
geno <- data.frame(geno, stringsAsFactors=FALSE, check.names=FALSE)
geno[geno=="0|0"] <- 0
geno[geno=="1|0" | geno=="0|1"] <- 1
geno[geno=="1|1"] <- 2
geno[geno==".|."] <- 3 # later replaced with NA 
# convert to numeric
geno[,10:ncol(geno)] <- sapply(geno[,10:ncol(geno)], as.numeric )
geno[geno==3] <- NA

maf <- unname(rowMeans(geno[,10:ncol(geno),with=FALSE])/2)
maf <- pmin(maf, 1-maf)

# transpose global covariates
ids <- globalcov[,ID]
globalcov[,ID := NULL]
subjs <- colnames(globalcov)
globalcov <- data.table(t(globalcov))
colnames(globalcov) <- ids
globalcov[,SUBJID:=subjs]
setkey(globalcov, SUBJID)

globalcov <- data.table(globalcov, keep.rownames="SUBJID")
setkey(globalcov, SUBJID)

allp <- file(out, 'w')
cat("gene_id\tvariant_id\ttss_distancet\tma_samples\tma_count\tmaf\tglobal_pval_nominal\tglobal_slope\tglobal_slope_se\tlava_pval_nominal\tlava_slope\tlava_slope_se\tafr_pval_nominal\tafr_slope\tafr_slope_se\n", file=allp)

# create a progress bar
total <- nrow(geno)*nrow(expr)
count = 0
system("echo 'Start eQTL calling with global and local covariates'")
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(j in 1:nrow(geno)){
	setTxtProgressBar(pb, count)
	ma_samples = length(grep('1|2', geno[j]))
	# make data table for local covariates
	variant_id = geno[j,ID]
	local <- localcov[SNP_ID == variant_id] # this should be the biggest time saver
	local <- local[,.(SUBJID,AFR,ASN,UNK)]
	strTmp = c('AFR','ASN','UNK') # Name of columns to be converted to numeric
	local[, (strTmp) := lapply(.SD, as.integer), .SDcols = strTmp]
	setkey(local, SUBJID)
	for(i in 1:nrow(expr)){
		count = count + 1
		gene_id = expr[i,gene_id]
		# check that the SNP is within 1 MB of the gene TSS 
		pos <- as.numeric(geno[j,POS])
		tss_distance = as.numeric(pos - expr[i,start])
		if(tss_distance <= 1000000 & tss_distance >= -1000000){
			# generate data frame for regression
			df1 <- data.table(expression=as.numeric(unname(unlist(expr[i]))[-c(1:4)]),haplo=as.numeric(unname(unlist(geno[j,]))[-c(1:9)]))
			df2 <- cbind(df1, globalcov)
			setkey(df2, SUBJID)
			dt <- df2[local,nomatch=0]
			strTmp = c('platform','sex','pcr')
			dt[, (strTmp) := lapply(.SD, as.factor), .SDcols = strTmp]
			ma_count = sum(dt[,haplo])

			df <- data.frame(dt, stringsAsFactors=FALSE)
			print(colnames(df))

			# regression with global PCs
			global.fit <- lm(expression ~ . - AFR - ASN - SUBJID, data=df)

			if(!is.null(warnings())){
				print(warnings())
			}
			out <- summary(global.fit)$coefficients

			if ("haplo"%in%rownames(out)==FALSE) {
				global_pval_nominal = "NA"
				global_slope = "NA"
				global_slope_se = "NA"
			} else {
				global_pval_nominal = out["haplo",4]
				global_slope = out["haplo",1]
				global_slope_se = out["haplo",2]
			}

			# LAVA
			local.fit <- lm(expression ~ . - PC1 - PC2 - PC3 - PC4 - PC5 - SUBJID, data=df)

			if(!is.null(warnings())){
				print(warnings())
			}
			out <- summary(local.fit)$coefficients

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
				lava_pval_nominal = "NA"
				lava_slope = "NA"
				lava_slope_se = "NA"
			} else {
				lava_pval_nominal = out["haplo",4]
				lava_slope = out["haplo",1]
				lava_slope_se = out["haplo",2]
			}
			m <- maf[j]

			string=paste(gene_id,variant_id,tss_distance,ma_samples,ma_count,m,global_pval_nominal,global_slope,global_slope_se,lava_pval_nominal,lava_slope,lava_slope_se,afr_pval_nominal,afr_slope,afr_slope_se, sep='\t')
			cat(paste(string,'\n',sep=''), file=allp, append=TRUE)
			break
		}
	}
}
close(pb)
close(allp)
