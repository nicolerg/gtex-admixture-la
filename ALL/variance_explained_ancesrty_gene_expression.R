# get variance explained by local ancestry versus global ancestry 

library(data.table)
library(variancePartition)

args <- commandArgs(trailingOnly=T)
tissue <- args[1]

expression = fread(cmd=sprintf('zcat /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/%s.v8.normalized_expression.bed.gz',tissue), sep='\t', header=T)
tss = fread(cmd='zcat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/hg19_tss_localcov.tsv.gz', sep='\t', header=T)
gpc = fread('/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_support_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt', sep='\t', header=T)
ax = fread('/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt',sep='\t',header=F)

gpc[,SUBJID := sapply(FID, function(x) paste(unname(unlist(strsplit(x, '-')))[1:2],collapse='-'))]
gpc = gpc[,.(SUBJID, PC1, PC2, PC3, PC4, PC5)]

# filter everything by admixed individuals with gene expression 
expression = expression[,c('gene_id',colnames(expression)[colnames(expression)%in%ax[,V1]]),with=F]
keep = colnames(expression)[2:ncol(expression)]
tissue
length(keep)
gpc = gpc[SUBJID %in% keep]
tss = tss[SUBJID %in% keep]

expression[,gene := sapply(gene_id, function(x) gsub("\\..*","",x))] # remove version 
expression[,gene_id := NULL]

# filter TSS by expressed genes
tss = tss[GENE_ID %in% expression[,gene]]

# remove rows with UNK
nrow(tss)
tss = tss[UNK == 0]
tss[,UNK := NULL]
nrow(tss)

merged_ancestry = merge(tss, gpc, by='SUBJID')

# fit a model for each gene and extract variance partition 

genes = expression[,gene]
expression[,gene:=NULL]
e = data.table(t(expression))
colnames(e) = genes 
e[,SUBJID := keep]

merged_ancestry[,ASN := as.numeric(ASN)]
merged_ancestry[,AFR := as.numeric(AFR)]

expr_var = file(sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/%s_var_explained_by_ancestry.tsv',tissue), 'w')
header = paste('gene', 'variance_explained_ga', 'variance_explained_la', sep='\t')
cat(paste0(header,'\n'), file=expr_var)

gene_ancestry_cor = file(sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/%s_var_ancestry_correlation.tsv',tissue), 'w')
header = paste('gene', 'PC1_AFR', 'PC1_ASN', 'PC2_AFR','PC2_ASN','PC3_AFR','PC3_ASN','PC4_AFR','PC4_ASN','PC5_AFR','PC5_ASN', sep='\t')
cat(paste0(header,'\n'), file=gene_ancestry_cor)

for (g in genes){
	
	# make data frame
	sub_expr = e[,c(g, 'SUBJID'),with=F]
	sub_meta = merged_ancestry[GENE_ID == g]
	if(nrow(sub_meta)==0){
		next
	}

	m = merge(sub_expr, sub_meta, by='SUBJID')
	colnames(m)[2] = 'gene'

	# # fit 
	# fit = lm(gene ~ ASN + AFR + PC1 + PC2 + PC3 + PC4 + PC5, data=m)

	# # get variance partition
	# varpart = extractVarPart(list(fit))
	# v = unname(unlist(varpart))
	# ga_var = sum(v[3:7])
	# la_var = sum(v[1:2])

	fit = lm(gene ~ ASN + AFR, m)
	m[,resid_la := resid(fit)]

	fit = lm(gene ~ PC1+PC2+PC3+PC4+PC5, m)
	m[,resid_ga := resid(fit)]
	
	# variance explained by GA 
	ga_var = summary(lm(resid_la ~ PC1+PC2+PC3+PC4+PC5, m))$r.squared

	# variance explained by LA 
	la_var = summary(lm(resid_ga ~ AFR + ASN, m))$r.squared

	# get correlation between each PC and LA for this gene 
	string = paste(g, cor(m[,PC1],m[,AFR]),
		cor(m[,PC1],m[,ASN]),
		cor(m[,PC2],m[,AFR]),
		cor(m[,PC2],m[,ASN]),
		cor(m[,PC3],m[,AFR]),
		cor(m[,PC3],m[,ASN]),
		cor(m[,PC4],m[,AFR]),
		cor(m[,PC4],m[,ASN]),
		cor(m[,PC5],m[,AFR]),
		cor(m[,PC5],m[,ASN]), sep='\t')
	cat(paste0(string,'\n'), file=gene_ancestry_cor, append=TRUE)

	# write out results 
	string = paste(g, ga_var, la_var, sep='\t')
	cat(paste0(string,'\n'), file=expr_var, append=TRUE)

}
close(expr_var)
close(gene_ancestry_cor)

