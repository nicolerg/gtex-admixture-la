# get number of admixed individuals in each tissue
library(data.table)
library(ggplot2)

#pheno <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
# pheno <- fread(pheno, sep='\t', header=TRUE)
# # aa <- pheno[RACE==2|RACE==1,SUBJID] # AA or Asian
# aa <- pheno[RACE==2,SUBJID] # AA or Asian

admixed <- fread('/mnt/lab_data/montgomery/nicolerg/rfmix-admixed/combined/gtex-admixed0.9.txt', sep='\t', header=FALSE)
aa <- admixed[,V1]

flist <- dir(path = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices", 
	pattern = "normalized_expression.bed.gz", 
	full.names = TRUE)

t <- c()
n <- c()
frac <- c()
total <- c()
for (f in flist){
	if (!grepl("gz.tbi",f,fixed=TRUE)){
		x <- fread(cmd=sprintf("zcat %s",f),sep='\t',header=TRUE,nrows=3)
		tissue <- gsub(".*/","",f)
		tissue <- gsub(".v8.normalized_expression.bed.gz","",tissue,fixed=TRUE)
		tot <- ncol(x) - 4
		print(tissue)
		subjs <- colnames(x)
		subjs <- subjs[subjs %in% aa]
		naa <- length(subjs)
		print(naa)
		t <- c(t, tissue)
		n <- c(n, naa)
		total <- c(total, tot)
		print(naa/tot)
		frac <- c(frac, round(naa/tot,2))
	}
}

print(min(frac))
print(max(frac))

df <- data.frame(tissue=t,N_AA=n,frac=frac)
df$frac <- as.character(df$frac)
df$frac <- gsub("0.",".",df$frac,fixed=TRUE)

print(sum(total))
print(sum(n))

write.table(df,'/mnt/lab_data/montgomery/nicolerg/local-eqtl/aa-per-tissue.tsv', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
