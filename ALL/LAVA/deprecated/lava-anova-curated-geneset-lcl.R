library(data.table)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]
geneset <- args[2]
out1 <- args[3]

tissue <- "Artery_Tibial"

geno <- paste0("/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/gtex.MAF01.phased.chr",chr,".vcf.gz")
expr <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/",tissue,".v8.normalized_expression.bed.gz"
globalcov <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/",tissue,".v8.covariates.txt"
localcov <- paste0("/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr",chr,".hg19.localcov.tsv.gz")
pheno <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
#curated <- paste0("/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/eqtl.curated.geneset.chr",chr,".tsv")
outdir <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/geuvadis-curated/"
#out1 <- paste0(outdir,"lymphocytes-local-allpairs-MAF0-chr",chr,".tsv")

system("echo 'Loading files'")
geno <- fread(input=sprintf("zcat %s", geno), sep='\t', header=TRUE)
expr <- fread(input=sprintf("zcat %s", expr), sep='\t', header=TRUE)
curated <- fread(geneset, sep='\t', header=TRUE)
globalcov <- fread(globalcov, sep='\t', header=TRUE)
localcov <- fread(input=sprintf("zcat %s", localcov), sep='\t', header=TRUE)
pheno <- fread(pheno, sep='\t', header=TRUE)
#curated <- fread(curated, sep='\t', header=TRUE)

#colnames(globalcov) <- gsub('.','-',colnames(globalcov),fixed=TRUE)
#colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

# select curated genes. reduce to gene id base for compatibility with different versions
genes <- curated[,pid]
genes <- unlist(lapply(genes, gsub, pattern="[.].*", replacement=""))
expr[,gene_id := unlist(lapply(expr[,gene_id], gsub, pattern="[.].*", replacement=""))]
expr <- expr[gene_id %in% genes]
print(paste0("N curated genes: ",length(genes)), quote=FALSE)
print(paste0("N curated genes in GTEx expression: ",nrow(expr)), quote=FALSE)

system("echo 'Sorting subjects in all tables'")
# filter pheno by expr
subjs <- colnames(expr)[5:ncol(expr)]
pheno <- pheno[SUBJID %in% subjs]
# only keep admixed individuals
nadx <- length(grep("1|2", pheno[,RACE]))
print(paste0("N Asian: ", length(grep("1", pheno[,RACE]))), quote=FALSE)
print(paste0("N black: ", length(grep("2", pheno[,RACE]))), quote=FALSE)
# filter expression file
subjs <- c(pheno[ RACE == '1' | RACE == '2', SUBJID ])
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

# transpose global covariates
ids <- globalcov[,ID]
globalcov[,ID := NULL]
subjs <- colnames(globalcov)
globalcov <- data.table(t(globalcov))
colnames(globalcov) <- ids
globalcov[,SUBJID:=subjs]
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
	# make data table for local covariates
	variant_id <- geno$ID[j]
	local <- localcov[SNP_ID == variant_id] # this should be the biggest time saver
	alocal <- local[,.(SUBJID,AFR,ASN,UNK)]
	strTmp = c('AFR','ASN','UNK') # Name of columns to be converted to numeric
	alocal[, (strTmp) := lapply(.SD, as.integer), .SDcols = strTmp]
	setkey(alocal, SUBJID)
	print(alocal)
	q()
	for(i in 1:nrow(expr)){
		gene_id = expr[i,gene_id]
		count = count + 1
		# check that the SNP is within 1 MB of the gene TSS 
		pos <- as.numeric(geno$POS[j])
		tss_distance = as.numeric(pos - expr[i,start])
		if(tss_distance <= 1000000 & tss_distance >= -1000000){
			# generate data frame for regression
			df1 <- data.table(expression=as.numeric(unname(unlist(expr[i]))[-c(1:4)]),
				haplo=as.numeric(unname(unlist(geno[j,]))[-c(1:9)]))
			df2 <- cbind(df1, globalcov)
			setkey(df2, SUBJID)
			dt <- df2[alocal,nomatch=0]
			strTmp = c('platform','sex','pcr')
			dt[, (strTmp) := lapply(.SD, as.factor), .SDcols = strTmp]

			df <- data.frame(dt, stringsAsFactors=FALSE)

			print(df)
			q()
			# replace NA genotype with colmean
			df$haplo[is.na(df$haplo)] <- mean(df$haplo, na.rm=TRUE)

			# regression without local covariates
			global.fit <- lm(expression ~ . - AFR - ASN - UNK - SUBJID, data=df)

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

			if ("AFR"%in%rownames(out)==FALSE) {
				afr_pval = "NA"
				afr_slope = "NA"
				afr_slope_se = "NA"
			} else {
				afr_pval = out["AFR",4]
				afr_slope = out["AFR",1]
				afr_slope_se = out["AFR",2]
			}

			if ("ASN"%in%rownames(out)==FALSE) {
				asn_pval = "NA"
				asn_slope = "NA"
				asn_slope_se = "NA"
			} else {
				asn_pval = out["ASN",4]
				asn_slope = out["ASN",1]
				asn_slope_se = out["ASN",2]
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
}
close(pb)
close(allp)
