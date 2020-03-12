library(data.table)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]
tissue <- args[2]
exprfile <- args[3]
out <- args[4]

geno <- paste0("/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/gtex.MAF01.phased.chr",chr,".vcf.gz")
globalcov <- paste0("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/",tissue,".v8.covariates.txt")
localcov <- paste0("/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr",chr,".hg19.localcov.tsv.gz")
pheno <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
globalai <- "/mnt/lab_data/montgomery/nicolerg/rfmix-outputs/gtex_ai_global.txt"

system("echo 'Loading files'")
geno <- fread(input=sprintf("zcat %s", geno), sep='\t', header=TRUE)
expr <- fread(input=sprintf("zcat %s", exprfile), sep='\t', header=TRUE)
globalcov <- fread(globalcov, sep='\t', header=TRUE)
localcov <- fread(input=sprintf("zcat %s", localcov), sep='\t', header=TRUE)
pheno <- fread(pheno, sep='\t', header=TRUE)

rows <- paste0("InferredCov",31:60)
globalcov <- globalcov[!ID %in% rows]

colnames(globalcov) <- gsub('.','-',colnames(globalcov),fixed=TRUE)
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

system("echo 'Sorting subjects in all tables'")
# filter pheno by expr
subjs <- colnames(expr)[5:ncol(expr)]
pheno <- pheno[SUBJID %in% subjs]
# filter expression file
subjs <- pheno[ RACE == '1' | RACE == '2', SUBJID ]

cols <- c("#chr","start","end","gene_id",subjs)
expr <- expr[,cols,with=FALSE]
# filter geno file
geno <- geno[,c(colnames(geno)[1:9], subjs), with=FALSE]
# filter global cov file
globalcov <- globalcov[,c("ID",subjs),with=FALSE]

system("echo 'Filtering variants'")
# remove variants that are too far from gene set
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

setDT(geno)

system("echo 'Formatting global covariates table'")
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

trylmglobal <- function(df) {
    out <- tryCatch(
        {
            # Just to highlight: if you want to use more than one 
            # R expression in the "try" part then you'll have to 
            # use curly brackets.
            # 'tryCatch()' will return the last evaluated expression 
            # in case the "try" part was completed successfully
			lm(expression ~ . - AFR - ASN - SUBJID, data=df)

            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            # message(cond)
            # Choose a return value in case of error
            return(NA)
        }
    )    
    return(out)
}

trylmlocal <- function(df) {
    out <- tryCatch(
        {
			lm(expression ~ . - PC1 - PC2 - PC3 - PC4 - PC5 - SUBJID, data=df)
        },
        error=function(cond) {
            # message(cond)
            return(NA)
        }
    )    
    return(out)
}

# create a progress bar
total <- nrow(geno)*nrow(expr)
count = 0
system("echo 'Start eQTL calling with global and local covariates'")
pb <- txtProgressBar(min = 0, max = total, style = 3)

for(j in 1:nrow(geno)){
	setTxtProgressBar(pb, count)

	# skip if no samples have this genotype
	haplo = as.numeric(unname(unlist(geno[j,]))[-c(1:9)])
	if (length(grep('1|2', haplo)) == 0){
		next
	}

	# make data table for local covariates
	variant_id = geno[j,ID]
	local <- localcov[SNP_ID == variant_id] # this should be the biggest time saver
	local <- local[ UNK == 0 ]
	local <- local[,.(SUBJID,AFR,ASN)]
	strTmp = c('AFR','ASN') # Name of columns to be converted to numeric
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
			df1 <- data.table(expression=as.numeric(unname(unlist(expr[i]))[-c(1:4)]),haplo=haplo)
			df2 <- cbind(df1, globalcov)
			setkey(df2, SUBJID)
			dt <- df2[local,nomatch=0]
			strTmp = c('platform','sex','pcr')
			dt[, (strTmp) := lapply(.SD, as.factor), .SDcols = strTmp]

			df <- data.frame(dt, stringsAsFactors=FALSE)
			df <- na.omit(df)

			ma_samples = length(grep('1|2', dt$haplo))
			ma_count = sum(dt$haplo)

			if(is.na(ma_count)){
				print(df)
			} else {
				next
			}

			if(length(unique(df$haplo))==1){
				next
			}

			# regression with global PCs
			global.fit <- trylmglobal(df)

			if(!is.na(global.fit[1])){
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
			} else {
				global_pval_nominal = "ERROR"
				global_slope = "ERROR"
				global_slope_se = "ERROR"
			}

			# LAVA
			local.fit <- trylmlocal(df)

			if(!is.na(local.fit[1])){
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
			} else {
				afr_pval_nominal = "ERROR"
				afr_slope = "ERROR"
				afr_slope_se = "ERROR"
				lava_pval_nominal = "ERROR"
				lava_slope = "ERROR"
				lava_slope_se = "ERROR"
			}

		}
		setTxtProgressBar(pb, count)
	}
}
