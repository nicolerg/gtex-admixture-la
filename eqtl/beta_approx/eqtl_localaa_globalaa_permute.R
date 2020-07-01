#!/bin/R

# this script assumes I've already run eQTL calling and pulls info from egenes_master
# this could easily be adapted to do both in one step 

library(data.table)
library(fitdistrplus)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1] # I don't actually use this in the script
tissue <- args[2]
exprfile <- args[3]
globalcov <- args[4]
geno <- args[5] # path to 117AX-filtered VCF for this chromosome
localcov <- args[6] # path to local ancestry covariate file for this chromosome
out <- args[7] # output file 
egenes_master_file <- args[8]
keep <- args[9]

timestamp = function(message){
	writeLines(sprintf('%s	%s', as.character(Sys.time()), message))
}

#              gene_id count              variant_id pval_nominal     slope
# 1: ENSG00000242268.2  4868  chr3_168698689_C_T_b38 0.0055884102 1.4030475
# 2: ENSG00000242268.2  4868  chr3_168699742_C_T_b38 0.0055884102 1.4030475
# 3: ENSG00000242268.2  4868  chr3_168706785_C_T_b38 0.0055884102 1.4030475
# 4: ENSG00000280143.1  6411 chr11_118266030_G_T_b38 0.0002753641 0.8618228
# 5: ENSG00000280143.1  6411 chr11_118266156_G_A_b38 0.0002753641 0.8618228
# 6: ENSG00000280143.1  6411 chr11_118266845_A_G_b38 0.0002753641 0.8618228
#    method               tissue
# 1: global Adipose_Subcutaneous
# 2: global Adipose_Subcutaneous
# 3: global Adipose_Subcutaneous
# 4: global Adipose_Subcutaneous
# 5: global Adipose_Subcutaneous
# 6: global Adipose_Subcutaneous

# subset down to this tissue 
load(egenes_master_file)
Tissue = tissue
egenes_master = egenes_master[tissue == Tissue]

timestamp('Loading files')
geno <- fread(cmd=sprintf("zcat %s", geno), sep='\t', header=TRUE) # filtered by chrom
expr <- fread(cmd=sprintf("zcat %s", exprfile), sep='\t', header=TRUE)
globalcov <- fread(globalcov, sep='\t', header=TRUE)
localcov <- fread(cmd=sprintf("zcat %s", localcov), sep='\t', header=TRUE)
keep <- fread(keep, sep='\t', header=FALSE)

N_ITER=100

rows <- paste0("InferredCov",31:60)
globalcov <- globalcov[!ID %in% rows]

colnames(globalcov) <- gsub('.','-',colnames(globalcov),fixed=TRUE)
colnames(geno) <- gsub('.','-',colnames(geno),fixed=TRUE)

timestamp('Sorting subjects in all tables')
# filter "keep" by expr
subjs <- colnames(expr)[5:ncol(expr)]
keep <- keep[V1 %in% subjs,V1]
# filter expression file
subjs <- subjs[na.omit(match(keep, subjs))]
cols <- c("#chr","start","end","gene_id",subjs)
expr <- expr[,cols,with=FALSE]
# filter geno file
geno <- geno[,c(colnames(geno)[1:9], subjs), with=FALSE]
setnames(geno, "ID", "VARIANT_ID")
# filter global cov file
globalcov <- globalcov[,c("ID",subjs),with=FALSE]

timestamp('Filtering variants')

# convert to hg38 position
extract_id <- function(x){
	return(as.integer(unname(unlist(strsplit(x, '_')))[2]))
}

geno[,POS := sapply(VARIANT_ID,extract_id)]

minstart <- min(as.numeric(expr[,start])) - 1000000 
maxstart <- max(as.numeric(expr[,start])) + 1000000
geno <- geno[POS >= minstart & POS <= maxstart]

localcov = localcov[SNP_ID %in% geno[,VARIANT_ID]]

timestamp('Converting phased haplotypes to genotype dosage')
convert_col = function(y){
	return(as.numeric(sapply(y, convert_hap)))
}
convert_hap = function(x){
	if (x == '0|0'){
		return('0')
	} else if (x == '0|1' | x == '1|0'){
		return('1')
	} else if (x == '1|1'){
		return('2')
	} else{
		return(NA)
	}
}
cols = colnames(geno)[grepl('GTEX', colnames(geno))]
geno[, (cols) := lapply(.SD, convert_col), .SDcols = cols]

timestamp('Formatting global covariates table')
# transpose global covariates
ids <- globalcov[,ID]
globalcov[,ID := NULL]
subjs <- colnames(globalcov)
globalcov <- data.table(t(globalcov))
colnames(globalcov) <- ids
globalcov[,SUBJID:=subjs]
setkey(globalcov, SUBJID)
strTmp = c('platform','sex','pcr')
globalcov[, (strTmp) := lapply(.SD, as.factor), .SDcols = strTmp]

trylmglobal <- function(df) {
    out = tryCatch(
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
    out = tryCatch(
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

timestamp('Permuting gene expression')
set.seed(1)
subjs = colnames(expr)[5:ncol(expr)]
permuted_expr = list()
permuted_pvalues = list()
for (gene in expr[,gene_id]){
	# permute gene expression 1000 times 
	permuted_expr[[gene]] = list()
	#permuted_pvalues[[gene]] = list()
	#permuted_pvalues[[gene]][['global']] = list()
	#permuted_pvalues[[gene]][['LAVA']] = list()
	e = as.numeric(unname(unlist(expr[gene_id==gene]))[5:ncol(expr)])
	for (i in 1:N_ITER){
		perm_e = sample(e, length(e), replace=F)
		df = data.frame(SUBJID = subjs, expression = perm_e)
		permuted_expr[[gene]][[i]] = df
		# permuted_pvalues[[gene]][['global']][[i]] = list()
		# permuted_pvalues[[gene]][['LAVA']][[i]] = list() # list of p-values for this permutation
	}
}

timestamp('Start tests for each SNP')

gc(verbose=F)


count = 0
total = nrow(geno)

allp = file(out, 'w')
header=paste('tissue','gene_id','permutation','variant_id','global_pval_nominal','local_pval_nominal', sep='\t')
cat(paste0(header,'\n'), file=allp)

for(j in 1:nrow(geno)){

	count = count + 1

	# skip if less than 5 samples have this genotype
	hap_dt = data.table(SUBJID = colnames(geno)[10:ncol(geno)], 
		haplo = as.numeric(unname(unlist(geno[j]))[-c(1:9)]))
	setkey(hap_dt, SUBJID)
	if(nrow(hap_dt[haplo > 0]) < 5){
		next
	}
	variant_id = geno[j, VARIANT_ID]
	timestamp(sprintf('Starting analysis for %s (%s/%s)',variant_id, count, total))

	# make data table for everything except expression 
	# do this only once per SNP 
	local = localcov[SNP_ID == variant_id] # this should be the biggest time saver
	local = local[ UNK == 0 ]
	local = local[,.(SUBJID,AFR,ASN)]
	strTmp = c('AFR','ASN') # Name of columns to be converted to numeric
	local[, (strTmp) := lapply(.SD, as.integer), .SDcols = strTmp]
	setkey(local, SUBJID)
	# global covariates 
	dt = globalcov[local,nomatch=0]
	setkey(dt, SUBJID)
	# genotype
	dt = dt[hap_dt, nomatch=0]
	
	df = data.frame(dt, stringsAsFactors=FALSE)
	df = na.omit(df)
	
	t2 = nrow(expr)
	c2 = 0
	pb = txtProgressBar(min = 0, max = t2, style = 3)

	for(i in 1:nrow(expr)){
	  c2 = c2 + 1
		gene_id = expr[i,gene_id]

		setTxtProgressBar(pb, c2)

		# check that the SNP is within 1 MB of the gene TSS 
		pos = as.numeric(geno[j,POS])
		tss_distance = as.numeric(pos - expr[i,start])
		
		if(tss_distance <= 1000000 & tss_distance >= -1000000){

			# test 1000 permutations of the gene expression 
			for (perm in 1:N_ITER){

				df2 = merge(df, permuted_expr[[gene_id]][[perm]], by='SUBJID')

				if(length(unique(df2$haplo))==1){
					next
				}

				# regression with global PCs
				global.fit = trylmglobal(df2)
				if(!is.na(global.fit[1])){
					global_pval = summary(global.fit)$coefficients["haplo",4]
				}else{
					global_pval = NA_real_
				}
				
				# LAVA
				local.fit <- trylmlocal(df2)
				if(!is.na(local.fit[1])){
					local_pval = summary(local.fit)$coefficients["haplo",4]
				}else{
					local_pval = NA_real_
				}
				string=paste(Tissue,gene_id,perm,variant_id,global_pval,local_pval, sep='\t')
				cat(paste0(string,'\n'), file=allp, append=T)

			}
		}
	}
	close(pb)
}
close(allp)

timestamp('Done generating null distributions of p-values')

outfile = sprintf('%s.pvalues.RData',gsub('\\..*', '',out))
save(permuted_expr, file=outfile)

# follow up with eqtl_localaa_globalaa_permute.R
