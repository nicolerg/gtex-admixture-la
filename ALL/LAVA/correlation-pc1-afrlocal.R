# change inputs to be standard outputs of LAVA 
allpairs <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/wholeblood-local-allpairs-chr21.tsv'

pc <- '/users/nicolerg/gtex-admix/correlation/pc1.tsv'
pvals <- '/users/nicolerg/gtex-admix/correlation/chr21-global-local.tsv'
afr <- '/users/nicolerg/gtex-admix/correlation/chr21-afr.tsv'

pc1 <- fread(input=pc, sep='\t', header=TRUE)
pval <- fread(input=pvals, sep='\t', header=TRUE)
geno <- fread(input=afr, sep='\t', header=TRUE)
pc1(setkey, ID)

unisnp <- unique(geno[SNP_ID])
correlations <- data.table(SNP=unisnp, PC1_AFRLOCAL_CORR=rep(NA, length(unisnp)), GLOBAL_LOCAL_PVAL_CORR=rep(NA, length(unisnp)))

counter = 0
paste('N variants:'length(unisnp),sep=' ')
for (snp in unisnp){
	counter = counter + 1
	print(paste('Working on SNP',snp,sep=' '))
	local <- geno[SNP_ID == snp, c(SUBJID,AFR)]
	setkey(local, SUBJID)
	dt <- local[pc1, match=0]
	localr2 = cor(dt[AFR],dt[PC1])
	ratios <- pval[variant_id == snp]
	pvalr2 <- cor(ratios[global_pval_nominal],ratios[local_pval_nominal])
	correlations[SNP == snp, PC1_AFRLOCAL_CORR] <- localr2
	correlations[SNP == snp, GLOBAL_LOCAL_PVAL_CORR] <- pvalr2
	print(paste('Completion (percent):',counter/length(unisnp)/100 sep=' '))
}

write.table(correlations, '/users/nicolerg/gtex-admix/correlation/pc1-afrdosage-correlations.tsv', sep='\t', col.names=TRUE, row.names=TRUE, quote=FALSE)