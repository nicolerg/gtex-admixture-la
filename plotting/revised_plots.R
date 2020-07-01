library(data.table)
library(ggplot2)
library(reshape)
library(OneR)
library(ggrepel)
library(argparse)
library(plyr)
library(ggcorrplot)
library(gridExtra)
library(gtable)
library(grid)

source('/oak/stanford/groups/smontgom/nicolerg/src/gtex-admixture-la/plotting/august_figures.R')

load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-1e-04-20200418.RData')
gene_names <- fread(sprintf('%s/gene_name_map.tsv',meta_dir), sep='\t', header=TRUE)
master_coloc[, gene_stable_id := sapply(gene_id, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)
master_coloc[,gwas := gwas_trait]
master_coloc[,gwas_trait := sapply(gwas, fix_label)]
mixed = c('Wojcik',
          'pgc.scz2',
          'RA_OKADA_TRANS_ETHNIC',
          'ISGC_Malik_2016_METASTROKE_all_strokes',
          'ILAE_Genetic_generalised_epilepsy',
          'DIAGRAM_T2D_TRANS_ETHNIC',
          'CARDIoGRAM_C4D_CAD_ADDITIVE')
master_coloc[,EUR := sapply(gwas_trait, function(x) ifelse(any(sapply(mixed, function(y) grepl(y, x))), 0, 1))]
master_coloc[,pheno := paste(tissue, gene_id, gwas_trait, sep=':')]

# ID subset of loci 
m_sub = copy(master_coloc)
m_sub[,coloc_sig := coloc_h4_global > 0.5 | coloc_h4_local > 0.5]
m_sub[,finemap_sig := finemap_clpp_global > 0.01 | finemap_clpp_local > 0.01]
nrow(m_sub[coloc_sig | finemap_sig])
m_sub = m_sub[coloc_sig & finemap_sig]

# but do they agree?
m_sub[,coloc_local_sig := as.numeric(coloc_h4_local > 0.5)]
m_sub[,coloc_global_sig := as.numeric(coloc_h4_global > 0.5)]
m_sub[,finemap_global_sig := as.numeric(finemap_clpp_global > 0.01)]
m_sub[,finemap_local_sig := as.numeric(finemap_clpp_local > 0.01)]

m_sub = m_sub[(coloc_local_sig*finemap_local_sig == 1) | (coloc_global_sig*finemap_global_sig == 1)]
nrow(m_sub)

m_sub[,coloc_global_better := ifelse(coloc_h4_global > coloc_h4_local, 1, 0)]
m_sub[,finemap_global_better := ifelse(finemap_clpp_global > finemap_clpp_local, 1, 0)]
table(m_sub[,finemap_global_better], m_sub[,coloc_global_better])

nrow(m_sub)
m_sub[,coloc_global_better := ifelse(coloc_h4_global > coloc_h4_local, 1, 0)]
m_sub[,finemap_global_better := ifelse(finemap_clpp_global > finemap_clpp_local, 1, 0)]
table(m_sub[finemap_sig==T,finemap_global_better])
m_sub = m_sub[(coloc_global_better == 1 & finemap_global_better == 1) | (coloc_global_better == 0 & finemap_global_better == 0)]
table(m_sub[coloc_sig==T,coloc_global_better])
# remove cases where both are significant 
m_sub = m_sub[!(coloc_h4_global > 0.5 & coloc_h4_local > 0.5)]
m_sub = m_sub[!(finemap_clpp_global > 0.01 & finemap_clpp_local > 0.01)]
nrow(m_sub) # 31

master_coloc[,unique_hit := ifelse(pheno %in% m_sub[,pheno], 1, 0)]

save(master_coloc, file='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-for-plotting.RData')
