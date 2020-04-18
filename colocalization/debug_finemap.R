library(data.table)
library(ggplot2)

setwd('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC')
load('master_coloc-1e-04-20200417.RData')

# pick sites with a COLOC hit 
master_coloc = master_coloc[coloc_h4_global > 0.5 | coloc_h4_local > 0.5]
master_coloc = unique(master_coloc[,.(gene_id, tissue, gwas_trait)])
nrow(master_coloc)
counts = data.table(table(master_coloc[,gwas_trait]))
counts = counts[order(N, decreasing = T)]

# ten jobs, for the top GWAS with counts 
gwas = counts[1:10,V1]
for(g in gwas){
  j = gsub('\\.hg38', '', g)
  overlap = list.files(path='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/overlap/coloc_tests',pattern = gsub('.','\\.',j,fixed=T),full.names = T)[1]
  overlap = fread(overlap, sep='\t', header=T)
  head(overlap)
  overlap[,tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(basename(x), '\\.filtered')))[1])]
  overlap[,gwas_trait := g]
  m = merge(master_coloc, overlap, by.x=c('gene_id','gwas_trait','tissue'),by.y=c('feature','gwas_trait','tissue'))
  write.table(m, file=sprintf('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/overlap/coloc_tests/test_1000G_vcf/%s.coloc-tests.txt',g), sep='\t', col.names=T, row.names=F, quote=F)
}
