#!/bin/R
# 21 May 2019

library(data.table)
library(ggplot2)

global_newpeer <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/new_peers/eqtl/Skin_Not_Sun_Exposed_Suprapubic/gtex.admix.global.egenes.tied.txt.gz'
local_newpeer <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/new_peers/eqtl/Skin_Not_Sun_Exposed_Suprapubic/gtex.admix.lava.egenes.tied.txt.gz'

global <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Skin_Not_Sun_Exposed_Suprapubic/gtex.admix.global.egenes.tied.txt.gz'
local <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Skin_Not_Sun_Exposed_Suprapubic/gtex.admix.lava.egenes.tied.txt.gz'

new_global <- fread(cmd=sprintf('zcat %s',global_newpeer), header=TRUE, sep='\t')
old_global <- fread(cmd=sprintf('zcat %s',global), header=TRUE, sep='\t')
new_global <- new_global[,.(gene_id, pval_nominal, slope)]
old_global <- old_global[,.(gene_id, pval_nominal, slope)]
global_merged <- merge(new_global, old_global, by='gene_id', suffixes=c('_new_peer','_old_peer'))

new_local <- fread(cmd=sprintf('zcat %s',local_newpeer), header=TRUE, sep='\t')
old_local <- fread(cmd=sprintf('zcat %s',local), header=TRUE, sep='\t')
new_local <- new_local[,.(gene_id, pval_nominal, slope)]
old_local <- old_local[,.(gene_id, pval_nominal, slope)]
local_merged <- merge(new_local, old_local, by='gene_id', suffixes=c('_new_peer','_old_peer'))

global_merged[,ancestry := 'global']
local_merged[,ancestry := 'local']

# get R^2 values
global_merged[,pval_nominal_old_peer := -log10(pval_nominal_old_peer)]
global_merged[,pval_nominal_new_peer := -log10(pval_nominal_new_peer)]
local_merged[,pval_nominal_old_peer := -log10(pval_nominal_old_peer)]
local_merged[,pval_nominal_new_peer := -log10(pval_nominal_new_peer)]

global_r2 <- summary(lm(pval_nominal_new_peer~pval_nominal_old_peer, global_merged))$r.squared
local_r2 <- summary(lm(pval_nominal_new_peer~pval_nominal_old_peer, local_merged))$r.squared

print(sprintf('Global R^2: %s',global_r2))
print(sprintf('Local R^2: %s',local_r2))

merged <- data.table(rbind(local_merged, global_merged))

pdf("~/gtex-admix/plots/final_figures/compare-peer-chr19.pdf", width=8, height=4)
ggplot(merged, aes(x=pval_nominal_old_peer, y=pval_nominal_new_peer)) +
	geom_point(alpha=0.7) +
	labs(x='GTEx v8 PEERs (-log10 p-value)', y='Admixed only PEERs (-log10 p-value)') +
	facet_wrap(~ ancestry, ncol = 2) +
	theme_bw()
dev.off()





