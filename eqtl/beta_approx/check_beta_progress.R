#!/bin/R

library(data.table)

# see how many genes on chr22 are near complete for beta approximation 
indir = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE'
setwd(indir)

test_num = fread(cmd='zcat chr22_eqtl_tests.gz', sep='\t', header=F)
colnames(test_num) = c('gene_id','n_snps','variant_id','tissue','method')
test_num[,variant_id := NULL]
test_num[,method := NULL] # always same number of SNPs tested in both genes 
test_num = unique(test_num)
# g = test_num[method=='global']
# l = test_num[method=='local']
# g[,method := NULL] 
# l[,method := NULL] 
# test_num = merge(g, l, by = c('gene_id','tissue'), suffixes = c('_global','_local'))
# all(test_num[,n_snps_global] == test_num[,n_snps_local])

# tissue	gene_id	permutation	variant_id	global_pval_nominal	local_pval_nominal
beta = fread('beta_tested_snps_new.tsv', sep='\t', header=F)
colnames(beta) = c('tissue','gene_id','variant_id')
head(beta)
beta = unique(beta)
dt = unique(beta[,.(tissue, gene_id)])
for (i in 1:nrow(dt)){
  print(i)
  snps = beta[gene_id == dt[i,gene_id] & tissue == dt[i, tissue], variant_id]
  dt[i, n_tested := length(unique(snps))]
}
head(dt)

# how many genes are done or close to done?
m = merge(test_num, dt, by=c('gene_id','tissue'))
m[,diff := n_snps - n_tested]
hist(m[,diff], breaks=500)
m[n_snps < 1000]
max(m[,diff]) # 8198

m[tissue == 'Muscle_Skeletal' & gene_id == 'ENSG00000228839.5']
m[diff < 0]

m = m[n_snps > 1000]
m = m[diff < 0.2*n_snps]
head(m)
nrow(m)
write.table(m, 'genes_to_include_in_beta.tsv', sep='\t', col.names=T, row.names=F, quote=F)
