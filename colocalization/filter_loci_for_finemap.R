library(data.table)
library(ggplot2)

setwd('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC')
load('master_coloc-1e-04-20200417.RData')

# pick sites I care about 
sig = master_coloc[coloc_h4_global > 0.5 | coloc_h4_local > 0.5] # coloc hit
sig = unique(sig[,.(gene_id, tissue, gwas_trait)])

insig = master_coloc[coloc_h4_global < 0.5 & coloc_h4_local < 0.5] # no coloc hit 
insig = unique(insig[,.(gene_id, tissue, gwas_trait)])

gwas = unique(master_coloc[,gwas_trait])
gwas = gwas[!is.na(gwas)]

# break up coloc hits versus non coloc hits. run coloc hits first 
for(g in gwas){
  j = gsub('\\.hg38', '', g)
  overlap = list.files(path='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/overlap/coloc_tests',pattern = gsub('.','\\.',j,fixed=T),full.names = T)[1]
  overlap = fread(overlap, sep='\t', header=T)
  overlap[,tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(basename(x), '\\.filtered')))[1])]
  overlap[,gwas_trait := g]
  
  # insig 
  insig_loci = merge(overlap, insig, by.y=c('gene_id','gwas_trait','tissue'),by.x=c('feature','gwas_trait','tissue'))
  write.table(insig_loci, file=sprintf('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/overlap/coloc_tests/rerun_1000G_vcf/insig/%s.coloc-tests.txt',g), sep='\t', col.names=T, row.names=F, quote=F)
  
  # sig 
  sig_loci = merge(overlap, sig, by.y=c('gene_id','gwas_trait','tissue'),by.x=c('feature','gwas_trait','tissue'))
  if(nrow(sig_loci)>0){
    write.table(sig_loci, file=sprintf('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/overlap/coloc_tests/rerun_1000G_vcf/sig/%s.coloc-tests.txt',g), sep='\t', col.names=T, row.names=F, quote=F)
  }
}
