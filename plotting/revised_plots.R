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
load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-1e-04-20200414.RData')

# master_coloc[,line_colour := ifelse(clpp_h4_global > clpp_h4_local, 'global', 'LAVA')]
#master_coloc[, label:= ifelse((clpp_h4_global < 0.60 & clpp_h4_local > 0.75) | (clpp_h4_local < 0.60 & clpp_h4_global > 0.75), paste0(feature,':',gwas_trait), NA)]

#master_coloc[,label := NULL]

# # remove genes that don't have different lead SNPs
# load("/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/mx_merged-all-annotations-v2.RData")
# la <- mx[,.(gene_id, tissue, lead_snp_D)]
# master_coloc <- merge(master_coloc, la, by.x=c('feature','tissue'), by.y=c('gene_id','tissue'))
# master_coloc <- master_coloc[lead_snp_D < 1]

gene_names <- fread(sprintf('%s/gene_name_map.tsv',meta_dir), sep='\t', header=TRUE)

master_coloc[, gene_stable_id := sapply(gene_id, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)

master_coloc[,gwas_trait := sapply(gwas_trait, fix_label)]
master_coloc[,label_global:= ifelse(clpp_h4_local < 0.50 & clpp_h4_global > 0.5 & (clpp_h4_global - clpp_h4_local) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
master_coloc[,label_local:= ifelse(clpp_h4_local > 0.5 & clpp_h4_global < 0.5 & (clpp_h4_local - clpp_h4_global) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

writeLines(unique(master_coloc[!is.na(label_local)|!is.na(label_global),gwas_trait]))
#master_coloc[,label := ifelse(((clpp_h4_local < 0.60 & clpp_h4_global > 0.75) | (clpp_h4_local > 0.75 & clpp_h4_global < 0.60)), gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]



master_coloc = master_coloc[grepl('Wojcik',gwas_trait) & tissue == 'Nerve_Tibial']

g <- ggplot() +
  
  geom_rect(data=data.frame(xmin=0.5,xmax=1,ymin=0,ymax=0.5), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=0.3, fill=unname(methodcols['global'])) + # global is better
  geom_rect(data=data.frame(xmin=0,xmax=0.5,ymin=0.5,ymax=1), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=0.3, fill=unname(methodcols['LAVA'])) + # local is better
  
  geom_rect(data=data.frame(xmin=0.5,xmax=1,ymin=0,ymax=0.5), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), colour=unname(methodcols['global']),fill=NA) + # global is better
  geom_rect(data=data.frame(xmin=0,xmax=0.5,ymin=0.5,ymax=1), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), colour=unname(methodcols['LAVA']),fill=NA) + # local is better
  
  geom_segment(data=data.frame(x=0,xend=1,y=0,yend=1),aes(x=x,y=y,xend=xend,yend=yend),colour='black',linetype='dashed') +
  geom_segment(data=data.frame(x=0,xend=1,y=0,yend=0),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
  geom_segment(data=data.frame(x=0,xend=1,y=1,yend=1),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
  geom_segment(data=data.frame(x=0,xend=0,y=1,yend=0),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
  geom_segment(data=data.frame(x=1,xend=1,y=0,yend=1),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
  
  geom_point(data=master_coloc, aes(x=clpp_h4_global, y=clpp_h4_local, colour=tissue),alpha=0.7) +
  geom_point(data=master_coloc[!is.na(label_local) | !is.na(label_global)],alpha=1,aes(fill=tissue,clpp_h4_global, y=clpp_h4_local),colour='black',shape=21) +
  
  # local better, fill
  geom_label_repel(data=master_coloc[!is.na(label_local)],
                   nudge_x=1.05 - master_coloc[!is.na(label_local),clpp_h4_global],
                   direction="y",
                   nudge_y=0.1,
                   hjust=0,
                   ylim=c(0.68,1),
                   aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local,fill=tissue,colour=tissue),
                   seed=1,
                   alpha=0.4,
                   size=3,
                   show.legend=F) + 
  
  # local better, no fill
  geom_label_repel(data=master_coloc[!is.na(label_local)],
                   nudge_x=1.05 - master_coloc[!is.na(label_local),clpp_h4_global],
                   direction="y",
                   nudge_y=0.1,
                   hjust=0,
                   ylim=c(0.68,1),
                   aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local,colour=tissue),
                   seed=1,
                   fill=NA,
                   size=3,
                   show.legend=F) + 
  
  # local better, black text 
  geom_label_repel(data=master_coloc[!is.na(label_local)],
                   nudge_x=1.05 - master_coloc[!is.na(label_local),clpp_h4_global],
                   direction="y",
                   nudge_y=0.1,
                   hjust=0,
                   ylim=c(0.68,1),
                   aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local),
                   seed=1,
                   fill=NA,
                   colour='black',
                   label.size=NA,
                   segment.size=NA,
                   size=3,
                   show.legend=F) + 
  
  geom_label_repel(data=master_coloc[!is.na(label_global)],
                   nudge_x=1.05 - master_coloc[!is.na(label_global),clpp_h4_global],
                   nudge_y=-0.3,
                   direction="y",
                   hjust=0,
                   aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global,fill=tissue,colour=tissue),
                   seed=1,
                   alpha=0.4,
                   ylim=c(0,0.68),
                   size=3,
                   show.legend=F) + 
  
  geom_label_repel(data=master_coloc[!is.na(label_global)],
                   nudge_x=1.05 - master_coloc[!is.na(label_global),clpp_h4_global],
                   nudge_y=-0.3,
                   direction="y",
                   hjust=0,
                   aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global,colour=tissue),
                   seed=1,
                   fill=NA,
                   ylim=c(0,0.68),
                   size=3,
                   show.legend=F) + 
  
  geom_label_repel(data=master_coloc[!is.na(label_global)],
                   nudge_x=1.05 - master_coloc[!is.na(label_global),clpp_h4_global],
                   nudge_y=-0.3,
                   direction="y",
                   hjust=0,
                   aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global),
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
  
  labs(x='GlobalAA colocalization probability (PP4)', y='LocalAA colocalization probability (PP4)') +
  theme_classic() +
  
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1.85), expand = c(0,0)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1), expand = c(0,0)) +
  guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=3),ncol=6),fill=FALSE,label=F) + 
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

png(paste0(plot_dir,'/fig4a_COLOC-local-global-p4-scatter-bytissue-labelled-PAGE-nerve.png'), width=11,height=8,unit='in',res=300)
g
dev.off()

master_coloc[clpp_h4_global > 0.75]

#################################################################################################



# fst regression

# locus_plots.R