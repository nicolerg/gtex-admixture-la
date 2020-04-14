library(data.table)
library(ggplot2)

setwd('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC')

load('master_coloc-1e-04-20200410.RData') # filtered down to loci of interest (eGenes at 1e-4; different lead SNPs)

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

master_coloc[,c('clpp_h0_global','clpp_h1_global','clpp_h2_global','clpp_h3_global',
                'clpp_h0_local','clpp_h1_local','clpp_h2_local','clpp_h3_local') := NULL]

g = finemap[method=='global']
l = finemap[method=='local']
l[,method:=NULL]
g[,method:=NULL]
f = merge(g, l, by=c('ref_snp','gene_id','tissue','gwas_trait'), suffixes = c('_global','_local'), all=T)

# in some cases, the same gene is tested for coloc twice
# take the coloc where the most SNPs were tested
# n_snps_global usually equals n_snps_local
m = f[,list(n_snps=max(n_snps_global),
            ref_snp = ref_snp[which.max(n_snps_global)],
            clpp_global = clpp_global[which.max(n_snps_global)],
            clpp_local = clpp_local[which.max(n_snps_global)],
            `-log_gwas_pval_global` = `-log_gwas_pval_global`[which.max(n_snps_global)],
            `-log_eqtl_pval_global` = `-log_eqtl_pval_global`[which.max(n_snps_global)],
            `-log_gwas_pval_local` = `-log_gwas_pval_local`[which.max(n_snps_global)],
            `-log_eqtl_pval_local` = `-log_eqtl_pval_local`[which.max(n_snps_global)],
            clpp_mod_global = clpp_mod_global[which.max(n_snps_global)],
            clpp_mod_local = clpp_mod_local[which.max(n_snps_global)]
            ), 
           by = c('tissue','gwas_trait','gene_id') ]

# are there duplicate tests for FINEMAP?
nrow(m) 
nrow(unique(m)) 
nrow(unique(m[,.(gene_id, tissue, gwas_trait)])) 

# merge with COLOC
setnames(master_coloc, "clpp_h4_global", "coloc_h4_global")
setnames(master_coloc, "clpp_h4_local", "coloc_h4_local")
setnames(m, "clpp_global", "finemap_clpp_global")
setnames(m, "clpp_local", "finemap_clpp_local")
setnames(m, "clpp_mod_global", "finemap_clpp_mod_global")
setnames(m, "clpp_mod_local", "finemap_clpp_mod_local")
all_merged = merge(master_coloc, m, by=c('gene_id','tissue','gwas_trait'), suffixes=c('_coloc','_finemap'), all.x=T)

library(ggrepel)

# look at correlation between COLOC and FINEMAP 
pdf('~/global_coloc_compare.pdf',width=12,height=8)
ggplot(all_merged, aes(x=coloc_h4_global, y=-log10(finemap_clpp_global))) +
  geom_point(alpha=0.3, aes(colour=tissue)) + 
  geom_abline(linetype='dashed') +
  theme_classic() +
  labs(x='COLOC pp4',y='FINEMAP CLPP (-log10)',title='GlobalAA') +
  geom_text_repel(data=all_merged[-log10(finemap_clpp_global) > 4.2],
                  aes(label=gwas_trait),
                  nudge_x=0.2,
                  hjust=0)
dev.off()

all_merged[-log10(finemap_clpp_global) > 4.2]

# wow that's awful 
# what are some of these loci that have high probabilities of colocalization with only one method?
master_coloc = all_merged
ggplot(master_coloc, aes(x=coloc_h4_global, y=coloc_h4_local, colour=tissue)) +
  
  geom_point() +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  
  labs(x='GlobalAA colocalization probability (PP4)', y='LocalAA colocalization probability (PP4)') +
  theme_classic() +
  #guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=3),ncol=6),fill=FALSE,label=F) + 
  theme(panel.grid=element_blank(),
        legend.position='top',
        legend.justification='left',
        legend.title=element_blank(),
        legend.margin=margin(b = -2, unit='mm'),
        legend.direction='horizontal',
        axis.title.x=element_text(hjust=0.15),
        legend.text=element_text(size=11)
        # axis.title=element_text(size=14)
  ) +
  annotate("text",x=0.75,y=0.025,hjust=0.5,label='Stronger GlobalAA coloc.') +
  annotate("text",x=0.25,y=0.975,hjust=0.5,label='Stronger LocalAA coloc.') 



