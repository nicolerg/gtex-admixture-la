library(data.table)
library(ggplot2)
library(ggrepel)

## DEPRECATED
# see format_coloc_v2.R

setwd('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC')

if(!file.exists('all_merged-egenes_coloc_finemap-1e-04-20200415.RData')){
  
  load('master_coloc-1e-04-20200414.RData') # filtered down to loci of interest (eGenes at 1e-4; different lead SNPs)
  
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
  # pick best per gene/tissue/trait combination 
  l = l[,list(ref_snp = ref_snp[which.min(clpp)],
              n_snps = n_snps[which.min(clpp)],
              clpp = min(clpp, na.rm=T),
              gwas_pval = `-log_gwas_pval`[which.min(clpp)],
              eqtl_pval = `-log_eqtl_pval`[which.min(clpp)]),
        by=c('gene_id','tissue','gwas_trait')]
  g = g[,list(ref_snp = ref_snp[which.min(clpp)],
              n_snps = n_snps[which.min(clpp)],
              clpp = min(clpp, na.rm=T),
              gwas_pval = `-log_gwas_pval`[which.min(clpp)],
              eqtl_pval = `-log_eqtl_pval`[which.min(clpp)]),
        by=c('gene_id','tissue','gwas_trait')]
  
  f = merge(g, l, by=c('gene_id','tissue','gwas_trait'), suffixes = c('_global','_local'), all=T)
  
  # are there duplicate tests for FINEMAP?
  nrow(f) 
  nrow(unique(f)) 
  nrow(unique(f[,.(gene_id, tissue, gwas_trait)])) 
  
  # merge with COLOC
  setnames(master_coloc, "clpp_h4_global", "coloc_h4_global")
  setnames(master_coloc, "clpp_h4_local", "coloc_h4_local")
  setnames(f, "clpp_global", "finemap_clpp_global")
  setnames(f, "clpp_local", "finemap_clpp_local")
  all_merged = merge(master_coloc, f, by=c('gene_id','tissue','gwas_trait'), suffixes=c('_coloc','_finemap'), all.x=T)
  table(all_merged[,ref_snp_local_finemap] == all_merged[,ref_snp_global_finemap])
  
  save(all_merged, file='all_merged-egenes_coloc_finemap-1e-04-20200415.RData')
}else{
  load('all_merged-egenes_coloc_finemap-1e-04-20200415.RData')
}

# we already know the correlation between COLOC and FINEMAP is going to be very low 
# # look at correlation between COLOC and FINEMAP 
# pdf('~/global_coloc_compare.pdf',width=12,height=8)
# ggplot(all_merged, aes(x=coloc_h4_global, y=-log10(finemap_clpp_global))) +
#   geom_point(alpha=0.3, aes(colour=tissue)) + 
#   geom_abline(linetype='dashed') +
#   theme_classic() +
#   labs(x='COLOC pp4',y='FINEMAP CLPP (-log10)',title='GlobalAA') +
#   geom_text_repel(data=all_merged[-log10(finemap_clpp_global) > 4.2],
#                   aes(label=gwas_trait),
#                   nudge_x=0.2,
#                   hjust=0)
# dev.off()
# 
# all_merged[-log10(finemap_clpp_global) > 4.2]
# 
# # wow that's awful 
# # what are some of these loci that have high probabilities of colocalization with only one method?
# master_coloc = all_merged
# ggplot(master_coloc, aes(x=coloc_h4_global, y=coloc_h4_local, colour=tissue)) +
#   
#   geom_point() +
#   geom_hline(yintercept = 0.5) +
#   geom_vline(xintercept = 0.5) +
#   
#   labs(x='GlobalAA colocalization probability (PP4)', y='LocalAA colocalization probability (PP4)') +
#   theme_classic() +
#   #guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=3),ncol=6),fill=FALSE,label=F) + 
#   theme(panel.grid=element_blank(),
#         legend.position='top',
#         legend.justification='left',
#         legend.title=element_blank(),
#         legend.margin=margin(b = -2, unit='mm'),
#         legend.direction='horizontal',
#         axis.title.x=element_text(hjust=0.15),
#         legend.text=element_text(size=11)
#         # axis.title=element_text(size=14)
#   ) +
#   annotate("text",x=0.75,y=0.025,hjust=0.5,label='Stronger GlobalAA coloc.') +
#   annotate("text",x=0.25,y=0.975,hjust=0.5,label='Stronger LocalAA coloc.') 

