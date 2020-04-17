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
  
  # split local and global 
  l = merged[method=='local']
  g = merged[method=='global']
  l = unique(l)
  g = unique(g)
  
  l[,method := NULL]
  g[,method := NULL]
  
  merged = merge(g, l, by=c('ref_snp','n_snps','tissue','gwas_trait','feature'), suffixes=c('_global','_local'), all=T)
  write.table(merged, 'merged_coloc_results.tsv', sep='\t', col.names=T, row.names=F)
}else{
  merged = fread('merged_coloc_results.tsv', sep='\t', header=T)
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
  
  g = finemap[method=='global']
  l = finemap[method=='local']
  l[,method:=NULL]
  g[,method:=NULL]
  
  merged_finemap = merge(g, l, by=c('ref_snp','n_snps','tissue','gwas_trait','gene_id'), suffixes=c('_global','_local'), all=T)
  write.table(merged_finemap, 'merged_finemap_results.tsv', sep='\t', col.names=T, row.names=F)
}else{
  merged_finemap = fread('merged_finemap_results.tsv', sep='\t', header=T)
}

# merge finemap and coloc
setnames(merged, "clpp_h4_global", "coloc_h4_global")
setnames(merged, "clpp_h4_local", "coloc_h4_local")
setnames(merged, "feature", "gene_id")
setnames(merged_finemap, "clpp_global", "finemap_clpp_global")
setnames(merged_finemap, "clpp_local", "finemap_clpp_local")
all_merged = merge(master_coloc, f, by=c('gene_id','tissue','gwas_trait'), suffixes=c('_coloc','_finemap'), all.x=T)
table(all_merged[,ref_snp_local_finemap] == all_merged[,ref_snp_global_finemap])

m = merge(merged, merged_finemap, by=c('ref_snp','tissue','gwas_trait','gene_id'), all=T, suffixes=c('_coloc','_finemap'))

# subset down to genes of interest
# merge with egenes 
load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/egenes_master-20200326.RData')
egenes = egenes_master[pval_nominal_global < 1e-4 & pval_nominal_local < 1e-4]
egenes = egenes[overlapping_lead_variants==0]
egenes = egenes[LD < 1]
egenes = egenes[,.(gene_id,tissue,pval_nominal_local,pval_nominal_global,lead_variants_global,lead_variants_local,LD)]
nrow(egenes)
master_coloc = merge(egenes,m, by=c('gene_id','tissue'), all.x=T)
nrow(master_coloc)

# when possible, keep the coloc result from the same seed snp 
indexes = (1:nrow(master_coloc))[!duplicated(master_coloc[,.(gene_id, tissue, gwas_trait)])]
keep_list = list()
x = 1 
for (i in indexes){
  sub = master_coloc[gene_id == gene_id[i] & tissue == tissue[i] & gwas_trait== gwas_trait[i]]
  sub = sub[ (is.na(`-log_eqtl_pval_global`) | `-log_eqtl_pval_global` > 4) & (is.na(`-log_eqtl_pval_local`) | `-log_eqtl_pval_local` > 4) ]
  if(nrow(sub) == 1){
    # just one line
    keep_list[[x]] = sub
  }else if(nrow(sub[complete.cases(sub)])>=1){
    # at least one case is complete 
    sub = sub[complete.cases(sub)]
    # pick ref_snp closest to global or local lead variant 
    # 7_8015011 chr7_7794838_A_G_b38
    sub[,pos_diff := mean(abs((as.numeric(unname(unlist(strsplit(ref_snp,'_')))[2]))-(as.numeric(unname(unlist(strsplit(unname(unlist(strsplit(lead_variants_local,';')))[1], '_')))[2]))), abs((as.numeric(unname(unlist(strsplit(ref_snp,'_')))[2]))-(as.numeric(unname(unlist(strsplit(unname(unlist(strsplit(lead_variants_global,';')))[1], '_')))[2]))))]
    sub = sub[which.min(pos_diff)]
    sub[,pos_diff := NULL]
    keep_list[[x]] = sub
  }else if(nrow(sub[complete.cases(sub)])==0){
    # no cases are complete 
    # first, narrow down to the most complete ref_snp(s)
    sub[,na_sum := is.na(coloc_h4_local) + is.na(coloc_h4_global) + is.na(finemap_clpp_global) + is.na(finemap_clpp_local)]
    sub = sub[which.min(na_sum)]
    sub[,na_sum := NULL]
    if(nrow(sub) > 1){
      sub[,pos_diff := mean(abs((as.numeric(unname(unlist(strsplit(ref_snp,'_')))[2]))-(as.numeric(unname(unlist(strsplit(unname(unlist(strsplit(lead_variants_local,';')))[1], '_')))[2]))), abs((as.numeric(unname(unlist(strsplit(ref_snp,'_')))[2]))-(as.numeric(unname(unlist(strsplit(unname(unlist(strsplit(lead_variants_global,';')))[1], '_')))[2]))))]
      sub = sub[which.min(pos_diff)]
      sub[,pos_diff := NULL]
    }
    keep_list[[x]] = sub
  }
  x = x + 1
}

pruned = rbindlist(keep_list)
pruned = pruned[,.(gene_id, tissue, ref_snp, gwas_trait, n_snps_coloc, coloc_h4_global, coloc_h4_local, 
                   n_snps_finemap, finemap_clpp_global, finemap_clpp_local, 
                   `-log_gwas_pval_global`, `-log_eqtl_pval_global`,
                   `-log_gwas_pval_local`, `-log_eqtl_pval_local`)]

p = merge(egenes, pruned, by=c('gene_id','tissue'), all.x=T)
nrow(unique(p[,.(gene_id, tissue, gwas_trait, ref_snp)])) - nrow(unique(p[,.(gene_id, tissue, gwas_trait)])) # 0 

master_coloc = p
head(master_coloc)
nrow(master_coloc)
nrow(master_coloc[complete.cases(master_coloc)])

# check for duplicates
master_coloc = unique(master_coloc)
nrow(master_coloc)
print(nrow(unique(master_coloc[,.(tissue,gwas_trait,gene_id)])))

save(master_coloc, file='master_coloc-1e-04-20200414.RData')
