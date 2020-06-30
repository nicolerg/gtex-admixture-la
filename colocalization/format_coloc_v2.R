#!/bin/R

# clean up coloc results 

library(data.table)

indir = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC'
setwd(indir)

if(!file.exists('merged_coloc_results.tsv')){
  
  fixed = fread(cmd=sprintf('zcat cc_fixed_coloc.tsv.gz'), sep='\t', header=T)
  nrow(fixed)
  head(fixed)
  # make tissue column
  fixed[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
  fixed[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
  table(fixed[,tissue])
  table(fixed[,method])
  fixed[,eqtl_file := NULL]
  fixed[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
  fixed[,c('base_gwas_file','gwas_trait'):=NULL]
  fixed = fixed[order(tissue, feature, gwas)]
  
  depr = fread(cmd=sprintf('zcat all_coloc_results.tsv.gz'), sep='\t', header=T)
  # format 
  depr[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
  depr[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
  depr[,eqtl_file := NULL]
  depr[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
  depr[,gwas := sapply(gwas, function(x) gsub('\\.txt\\.gz','',x))]
  depr[,gwas := sapply(gwas, function(x) gsub('coloc_imputed_','imputed_',x))]
  depr[,c('base_gwas_file','gwas_trait'):=NULL]
  depr = depr[order(tissue, feature, gwas)]
  
  # now remove gwas present in "fixed"
  gwas_fixed = unique(fixed[,gwas])
  length(unique(depr[,gwas]))
  depr = depr[!gwas %in% gwas_fixed]
  merged = rbind(depr, fixed)
  
  merged = merged[order(tissue, feature, gwas)]
  table(merged[,gwas])
  table(merged[,method], merged[,tissue])
  
  head(merged)
  setnames(merged, "gwas", "gwas_trait")
  
  # remove tests with low numbers of snps
  print(nrow(merged[n_snps <= 20]))
  merged = merged[n_snps > 20]
  
  # collapse by best coloc per gene/tissue/gwas/method combination 
  coloc_best = merged[,list(coloc_h4 = max(clpp_h4, na.rm=T),
                            ref_snp = ref_snp[which.max(clpp_h4)],
                            n_snps = n_snps[which.max(clpp_h4)]),
                      by=c('tissue','feature','gwas_trait','method')]
  coloc_best[coloc_h4 == -Inf, coloc_h4 := NA]
  
  # split local and global 
  l = coloc_best[method=='local']
  g = coloc_best[method=='global']
  l = unique(l[complete.cases(l)])
  g = unique(g[complete.cases(g)])
  
  l[,method := NULL]
  g[,method := NULL]

  merged = merge(g, l, by=c('tissue','gwas_trait','feature'), suffixes=c('_global','_local'), all=T)
  write.table(merged, 'merged_coloc_results.tsv', sep='\t', col.names=T, row.names=F, quote=F)
  merged_coloc = merged
}else{
  merged_coloc = fread('merged_coloc_results.tsv', sep='\t', header=T)
}

if(!file.exists('merged_finemap_results.tsv')){
  # now add in finemap
  finemap=fread(cmd='zcat all_finemap_results.tsv.gz', sep='\t', header=T)
  finemap[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
  finemap[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
  table(finemap[,tissue])
  table(finemap[,method])
  finemap[,eqtl_file := NULL]
  finemap[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
  length(unique(finemap[,gwas]))
  finemap[,c('base_gwas_file','gwas_trait'):=NULL]
  finemap = finemap[order(tissue, feature, gwas)]
  setnames(finemap, "gwas", "gwas_trait")
  setnames(finemap, "feature", "gene_id")
  
  # remove tests with low numbers of snps
  nrow(finemap)
  nrow(finemap[n_snps <= 20])
  finemap = finemap[n_snps > 20]
  
  # collapse by best coloc per gene/tissue/gwas/method combination 
  finemap_best = finemap[,list(finemap_clpp = max(clpp, na.rm=T),
                               finemap_clpp_mod = clpp_mod[which.max(clpp)],
                               gwas_pval_log10 = `-log_gwas_pval`[which.max(clpp)],
                               eqtl_pval_log10 = `-log_eqtl_pval`[which.max(clpp)],
                               ref_snp = ref_snp[which.max(clpp)],
                               n_snps = n_snps[which.max(clpp)]),
                         by=c('tissue','gene_id','gwas_trait','method')]
  
  g = finemap_best[method=='global']
  l = finemap_best[method=='local']
  l[,method:=NULL]
  g[,method:=NULL]
  l = unique(l[complete.cases(l)])
  g = unique(g[complete.cases(g)])
  
  merged_finemap = merge(g, l, by=c('tissue','gwas_trait','gene_id'), suffixes=c('_global','_local'), all=T)
  write.table(merged_finemap, 'merged_finemap_results.tsv', sep='\t', col.names=T, row.names=F, quote=F)
}else{
  merged_finemap = fread('merged_finemap_results.tsv', sep='\t', header=T)
}

# FINEMAP v2: 1000G VCF
finemap = fread(cmd='zcat signif_coloc_loci_tests_finemap_clpp_status.txt.gz', sep='\t', header=T)
finemap[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
finemap[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
table(finemap[,tissue])
table(finemap[,method])
finemap[,eqtl_file := NULL]
finemap[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
length(unique(finemap[,gwas]))
finemap[,c('base_gwas_file','gwas_trait'):=NULL]
finemap = finemap[order(tissue, feature, gwas)]
setnames(finemap, "gwas", "gwas_trait")
setnames(finemap, "feature", "gene_id")

# remove tests with low numbers of snps
nrow(finemap)
nrow(finemap[n_snps <= 20])
finemap = finemap[n_snps > 20]

# collapse by best coloc per gene/tissue/gwas/method combination 
finemap_best = finemap[,list(finemap_clpp = max(clpp, na.rm=T),
                             finemap_clpp_mod = clpp_mod[which.max(clpp)],
                             gwas_pval_log10 = `-log_gwas_pval`[which.max(clpp)],
                             eqtl_pval_log10 = `-log_eqtl_pval`[which.max(clpp)],
                             ref_snp = ref_snp[which.max(clpp)],
                             n_snps = n_snps[which.max(clpp)]),
                       by=c('tissue','gene_id','gwas_trait','method')]

g = finemap_best[method=='global']
l = finemap_best[method=='local']
l[,method:=NULL]
g[,method:=NULL]
l = unique(l[complete.cases(l)])
g = unique(g[complete.cases(g)])

merged_finemap = merge(g, l, by=c('tissue','gwas_trait','gene_id'), suffixes=c('_global','_local'), all=T)
write.table(merged_finemap, 'merged_finemap_results_1000G.tsv', sep='\t', col.names=T, row.names=F, quote=F)

# merge finemap and coloc
setnames(merged_coloc, "feature", "gene_id")
master_coloc = merge(merged_coloc, merged_finemap, by=c('gene_id','tissue','gwas_trait'), all=T, suffixes=c('_coloc','_finemap'))
write.table(master_coloc, file='merged_coloc_finemap_1000G_all.tsv', sep='\t', col.names=T, row.names=F, quote=F)

# subset down to genes of interest
# merge with egenes 
load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/egenes_master-20200326.RData')
egenes = egenes_master[pval_nominal_global < 1e-4 & pval_nominal_local < 1e-4]
egenes = egenes[overlapping_lead_variants==0]
egenes = egenes[LD < 1]
egenes = egenes[,.(gene_id,tissue,pval_nominal_local,pval_nominal_global,lead_variants_global,lead_variants_local,LD)]
nrow(egenes)
master_coloc = merge(egenes,master_coloc, by=c('gene_id','tissue'), all.x=T)
nrow(master_coloc)
nrow(unique(master_coloc[,.(gene_id,tissue)])) == nrow(egenes)

save(master_coloc, file='master_coloc-1e-04-20200418.RData')

hist(master_coloc[,finemap_clpp_global], breaks=500)

#####################################################################

# compare finemap GTEx vs 1KG VCF:
kg = fread(cmd='zcat signif_coloc_loci_tests_finemap_clpp_status.txt.gz',sep='\t',header=T)
kg[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
kg[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
table(kg[,tissue])
table(kg[,method])
kg[,eqtl_file := NULL]
kg[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
kg[,c('base_gwas_file','gwas_trait'):=NULL]
setnames(kg, "gwas", "gwas_trait")
setnames(kg, "feature", "gene_id")

# now add in finemap
finemap=fread(cmd='zcat all_finemap_results.tsv.gz', sep='\t', header=T)
finemap[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
finemap[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
table(finemap[,tissue])
table(finemap[,method])
finemap[,eqtl_file := NULL]
finemap[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
length(unique(finemap[,gwas]))
finemap[,c('base_gwas_file','gwas_trait'):=NULL]
finemap = finemap[order(tissue, feature, gwas)]
setnames(finemap, "gwas", "gwas_trait")
setnames(finemap, "feature", "gene_id")

f = merge(finemap, kg, by=c('tissue','gwas_trait','gene_id','ref_snp','method'), suffixes=c('_gtex','_1000G'))

library(ggplot2)

ggplot(f, aes(x=clpp_gtex, y=clpp_1000G)) +
  geom_point() +
  theme_classic() + 
  geom_abline(linetype='dashed')
