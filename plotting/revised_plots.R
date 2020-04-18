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

# coloc ################################################################################################

# already excludes loci with lead variants in LD of 1
#load(sprintf('%s/master_coloc-1e-04-20200414.RData',master_data_dir))
load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-1e-04-20200417.RData')

gene_names <- fread(sprintf('%s/gene_name_map.tsv',meta_dir), sep='\t', header=TRUE)

master_coloc[, gene_stable_id := sapply(gene_id, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)

master_coloc[,gwas_trait := sapply(gwas_trait, fix_label)]
master_coloc[,label_global:= ifelse(coloc_h4_local < 0.50 & coloc_h4_global > 0.5 & (coloc_h4_global - coloc_h4_local) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
master_coloc[,label_local:= ifelse(coloc_h4_local > 0.5 & coloc_h4_global < 0.5 & (coloc_h4_local - coloc_h4_global) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

df = unique(master_coloc[!is.na(label_local)|!is.na(label_global)])
write.table(df, '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/plots/coloc_diff.csv',sep=',',col.names=T,row.names=F,quote=F)

g <- ggplot() +

  geom_segment(data=data.frame(x=0,xend=1,y=0.5,yend=0.5),aes(x=x,y=y,xend=xend,yend=yend),colour=methodcols[['LAVA']]) +
  geom_segment(data=data.frame(x=0.5,xend=0.5,y=0,yend=1),aes(x=x,y=y,xend=xend,yend=yend),colour=methodcols[['global']]) +
  
  geom_point(data=master_coloc[coloc_h4_global < 0.5 & coloc_h4_local < 0.5], aes(x=coloc_h4_global, y=coloc_h4_local), colour='gray', alpha=0.3) +
  geom_point(data=master_coloc[coloc_h4_global > 0.5 | coloc_h4_local > 0.5], aes(x=coloc_h4_global, y=coloc_h4_local, colour=tissue),alpha=0.7) +
  geom_point(data=master_coloc[!is.na(label_local) | !is.na(label_global)],alpha=1,aes(fill=tissue, x=coloc_h4_global, y=coloc_h4_local),colour='black',shape=21) +
  
  # local better, fill
  geom_label_repel(data=master_coloc[!is.na(label_local)],
                   nudge_x=1.05 - master_coloc[!is.na(label_local),coloc_h4_global],
                   direction="y",
                   nudge_y=0.1,
                   hjust=0,
                   ylim=c(0.68,1),
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=label_local,fill=tissue,colour=tissue),
                   seed=1,
                   alpha=0.4,
                   size=3,
                   show.legend=F) + 
  
  # local better, no fill
  geom_label_repel(data=master_coloc[!is.na(label_local)],
                   nudge_x=1.05 - master_coloc[!is.na(label_local),coloc_h4_global],
                   direction="y",
                   nudge_y=0.1,
                   hjust=0,
                   ylim=c(0.68,1),
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=label_local,colour=tissue),
                   seed=1,
                   fill=NA,
                   size=3,
                   show.legend=F) + 
  
  # local better, black text 
  geom_label_repel(data=master_coloc[!is.na(label_local)],
                   nudge_x=1.05 - master_coloc[!is.na(label_local),coloc_h4_global],
                   direction="y",
                   nudge_y=0.1,
                   hjust=0,
                   ylim=c(0.68,1),
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=label_local),
                   seed=1,
                   fill=NA,
                   colour='black',
                   label.size=NA,
                   segment.size=NA,
                   size=3,
                   show.legend=F) + 
  
  geom_label_repel(data=master_coloc[!is.na(label_global)],
                   nudge_x=1.05 - master_coloc[!is.na(label_global),coloc_h4_global],
                   nudge_y=-0.3,
                   direction="y",
                   hjust=0,
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=label_global,fill=tissue,colour=tissue),
                   seed=1,
                   alpha=0.4,
                   ylim=c(0,0.68),
                   size=3,
                   show.legend=F) + 
  
  geom_label_repel(data=master_coloc[!is.na(label_global)],
                   nudge_x=1.05 - master_coloc[!is.na(label_global),coloc_h4_global],
                   nudge_y=-0.3,
                   direction="y",
                   hjust=0,
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=label_global,colour=tissue),
                   seed=1,
                   fill=NA,
                   ylim=c(0,0.68),
                   size=3,
                   show.legend=F) + 
  
  geom_label_repel(data=master_coloc[!is.na(label_global)],
                   nudge_x=1.05 - master_coloc[!is.na(label_global),coloc_h4_global],
                   nudge_y=-0.3,
                   direction="y",
                   hjust=0,
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=label_global),
                   seed=1,
                   fill=NA,
                   colour='black',
                   label.size=NA,
                   segment.size=NA,
                   ylim=c(0,0.68),
                   size=3,
                   show.legend=F) +
  
  scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue', limits=tissues) +
  scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +
  
  labs(x='GlobalAA COLOC PP4', y='LocalAA COLOC PP4') +
  theme_classic() +
  
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1.85), expand = c(0,0)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1), expand = c(0,0)) +
  guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=3),ncol=7),fill=FALSE,label=F) + 
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
  geom_abline(linetype='dashed')

png(paste0(plot_dir,'/fig4_COLOC.png'), width=9,height=6,unit='in',res=300)
g
dev.off()

#################################################################################################
load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-1e-04-20200417.RData')
gene_names <- fread(sprintf('%s/gene_name_map.tsv',meta_dir), sep='\t', header=TRUE)
master_coloc[, gene_stable_id := sapply(gene_id, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)
master_coloc[,gwas_trait := sapply(gwas_trait, fix_label)]

master_coloc[finemap_clpp_local > 0.01 | finemap_clpp_global > 0.01]
master_coloc = master_coloc[!is.na(finemap_clpp_local) & !is.na(finemap_clpp_global)]

master_coloc[order(finemap_clpp_mod_global, decreasing = T)][1:50]


max(master_coloc[,finemap_clpp_global])
max(master_coloc[,finemap_clpp_local])

g <- ggplot() +
  
  geom_segment(data=data.frame(x=0,xend=0.05,y=0.01,yend=0.01),aes(x=x,y=y,xend=xend,yend=yend),colour=methodcols[['LAVA']]) +
  geom_segment(data=data.frame(x=0.01,xend=0.01,y=0,yend=0.05),aes(x=x,y=y,xend=xend,yend=yend),colour=methodcols[['global']]) +
  
  geom_point(data=master_coloc[finemap_clpp_global<0.01 & finemap_clpp_local<0.01], aes(x=finemap_clpp_global, y=finemap_clpp_local), colour='gray', alpha=0.3) + 
  geom_point(data=master_coloc[!(finemap_clpp_global<0.01 & finemap_clpp_local<0.01)], aes(x=finemap_clpp_global, y=finemap_clpp_local, colour=tissue), alpha=0.7) + 

  scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue', limits=tissues) +
  scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +
  
  labs(x='GlobalAA FINEMAP CLPP', y='LocalAA FINEMAP CLPP') +
  theme_classic() +

  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=3),ncol=7),fill=FALSE,label=F) + 
  theme(panel.grid=element_blank(),
        legend.position='none',
        legend.justification='center',
        legend.title=element_blank(),
        legend.margin=margin(b = -2, unit='mm'),
        legend.direction='horizontal'
  ) +
  geom_abline(linetype='dashed') 

png(paste0(plot_dir,'/Fig4-FINEMAP.png'), width=4.5,height=3.5,unit='in',res=300)
g
dev.off()

#################################################################################################


# fst regression

# locus_plots.R