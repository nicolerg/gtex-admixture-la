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
master_coloc[,gwas_trait := sapply(gwas_trait, fix_label)]
mixed = c('Wojcik',
          'pgc.scz2',
          'RA_OKADA_TRANS_ETHNIC',
          'ISGC_Malik_2016_METASTROKE_all_strokes',
          'ILAE_Genetic_generalised_epilepsy',
          'DIAGRAM_T2D_TRANS_ETHNIC',
          'CARDIoGRAM_C4D_CAD_ADDITIVE')
master_coloc[,EUR := sapply(gwas_trait, function(x) ifelse(any(sapply(mixed, function(y) grepl(y, x))), 0, 1))]
master_coloc[,pheno := paste(tissue, gene_id, gwas_trait, sep=':')]

finemap_off_axis = master_coloc[((finemap_clpp_local < 0.01 & finemap_clpp_global > 0.01) | (finemap_clpp_local > 0.01 & finemap_clpp_global < 0.01)) & abs(finemap_clpp_local-finemap_clpp_global)>0.02, pheno ]
coloc_off_axis = master_coloc[((coloc_h4_local < 0.50 & coloc_h4_global > 0.5) | (coloc_h4_local > 0.50 & coloc_h4_global < 0.5))&abs(coloc_h4_global-coloc_h4_local) > 0.1, pheno ]
off_axis_both = intersect(finemap_off_axis, coloc_off_axis)

master_coloc[, label := ifelse(pheno %in% off_axis_both, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
master_coloc[!is.na(label), shortl := c('A','B','C','D','E','F')]

# coloc ################################################################################################

# already excludes loci with lead variants in LD of 1
#load(sprintf('%s/master_coloc-1e-04-20200414.RData',master_data_dir))

# master_coloc[,label_global:= ifelse(coloc_h4_local < 0.50 & coloc_h4_global > 0.5 & (coloc_h4_global - coloc_h4_local) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
# master_coloc[,label_local:= ifelse(coloc_h4_local > 0.5 & coloc_h4_global < 0.5 & (coloc_h4_local - coloc_h4_global) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
# df = unique(master_coloc[!is.na(label_local)|!is.na(label_global)])
#write.table(df, '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/plots/coloc_diff.csv',sep=',',col.names=T,row.names=F,quote=F)

g <- ggplot() +

  geom_hline(yintercept=0.5, colour=methodcols[['LAVA']]) +
  geom_vline(xintercept=0.5, colour=methodcols[['global']]) +
  
  geom_point(data=master_coloc[coloc_h4_global < 0.5 & coloc_h4_local < 0.5], aes(x=coloc_h4_global, y=coloc_h4_local, shape=factor(EUR)), fill='gray', colour='gray', alpha=0.3,size=1) +
  geom_point(data=master_coloc[coloc_h4_global > 0.5 | coloc_h4_local > 0.5], aes(x=coloc_h4_global, y=coloc_h4_local, fill=tissue, shape=factor(EUR), colour=tissue),alpha=0.7,size=1) +
  geom_point(data=master_coloc[!is.na(label)],alpha=1,aes(fill=tissue, x=coloc_h4_global, y=coloc_h4_local, shape=factor(EUR)),colour='black',size=2) +
  
  scale_shape_manual(values=c('1'=21,'0'=24)) +
  
  # fill
  geom_label_repel(data=master_coloc[!is.na(shortl)],
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=shortl,fill=tissue,colour=tissue),
                   nudge_x = 0.01,
                   segment.size=NA,
                   seed=1,
                   alpha=0.4,
                   size=3,
                   show.legend=F) + 
  
  # no fill
  geom_label_repel(data=master_coloc[!is.na(shortl)],
                   nudge_x = 0.01,
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=shortl,colour=tissue),
                   seed=1,
                   fill=NA,
                   segment.size=NA,
                   size=3,
                   show.legend=F) + 
  
  # black text 
  geom_label_repel(data=master_coloc[!is.na(shortl)],
                   nudge_x = 0.01,
                   aes(x=coloc_h4_global,y=coloc_h4_local,label=shortl),
                   seed=1,
                   fill=NA,
                   colour='black',
                   label.size=NA,
                   segment.size=NA,
                   size=3,
                   show.legend=F) + 
  
  scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue', limits=tissues) +
  scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +
  
  labs(x='GlobalAA COLOC PP4', y='LocalAA COLOC PP4') +
  theme_classic() +
  
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1), expand = c(0,0), labels=c('0','0.25','0.5','0.75','1')) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1), expand = c(0,0), labels=c('0','0.25','0.5','0.75','1')) +
  guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=2),ncol=4),
         fill=FALSE,
         label=F,
         shape = F) + 
  theme(panel.grid=element_blank(),
        legend.position = 'none',
        #legend.key.height = unit(1.5, 'lines'),
        #legend.justification='left',
        legend.title=element_blank(),
        legend.margin=margin(b = -2, unit='mm'),
        legend.direction='horizontal',
        legend.text=element_text(size=7)
        # axis.title=element_text(size=14)
  ) +
  geom_abline(linetype='dashed')

png(paste0(plot_dir,'/fig4_COLOC.png'), width=5,height=5,unit='in',res=300)
g
dev.off()

pdf(paste0(plot_dir,'/legend_only.pdf'), width=12,height=5)
ggplot() +
  #geom_point(data=master_coloc[coloc_h4_global < 0.5 & coloc_h4_local < 0.5], aes(x=coloc_h4_global, y=coloc_h4_local, shape=factor(EUR)), fill='gray', colour='gray', alpha=0.3,size=1) +
  geom_point(data=master_coloc[coloc_h4_global > 0.5 | coloc_h4_local > 0.5], aes(x=coloc_h4_global, y=coloc_h4_local, fill=tissue, shape=factor(EUR), colour=tissue),alpha=0.7,size=1) +
  geom_point(data=master_coloc[!is.na(label)],alpha=1,aes(fill=tissue, x=coloc_h4_global, y=coloc_h4_local, shape=factor(EUR)),colour='black',size=2) +
  
  scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue', limits=tissues) +
  scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +
  scale_shape_manual(values=c('1'=21,'0'=24), labels=c('1'='European','0'='Mixed ancestry')) +
  guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=2),ncol=7),
         shape = guide_legend(override.aes=list(alpha=1,size=2,colour='black',fill='black')),
         fill=FALSE,label=FALSE) + 
  theme_classic() +
  theme(panel.grid=element_blank(),
        legend.position = 'top',
        legend.title=element_blank(),
        legend.direction='horizontal',
        legend.text=element_text(size=7)) 
dev.off()

#################################################################################################
nrow(master_coloc[finemap_clpp_local > 0.01 | finemap_clpp_global > 0.01])
finemap = master_coloc[!is.na(finemap_clpp_local) & !is.na(finemap_clpp_global)]
max(finemap[,finemap_clpp_global])
max(finemap[,finemap_clpp_local])

#finemap[,label:= ifelse( ((finemap_clpp_local > 0.01 & finemap_clpp_global < 0.01) | (finemap_clpp_local < 0.01 & finemap_clpp_global > 0.03)) & abs(finemap_clpp_local - finemap_clpp_global) > 0.01, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

finemap[finemap_clpp_global < .0008, finemap_clpp_global := 0]
finemap[finemap_clpp_local < .0008, finemap_clpp_global := 0]

g2 <- ggplot() +
  
  geom_hline(yintercept=0.01, colour=methodcols[['LAVA']]) +
  geom_vline(xintercept=0.01, colour=methodcols[['global']]) +
  
  geom_point(data=finemap[finemap_clpp_global < 0.01 & finemap_clpp_local < 0.01], aes(x=finemap_clpp_global, y=finemap_clpp_local, shape=factor(EUR)), fill='gray', colour='gray', alpha=0.3,size=1) +
  geom_point(data=finemap[finemap_clpp_global > 0.01 | finemap_clpp_local > 0.01], aes(x=finemap_clpp_global, y=finemap_clpp_local, fill=tissue, shape=factor(EUR), colour=tissue),alpha=0.7,size=1) +
  geom_point(data=finemap[!is.na(label)],alpha=1,aes(fill=tissue, x=finemap_clpp_global, y=finemap_clpp_local, shape=factor(EUR)),colour='black',size=2) +
  
  scale_colour_manual(values=tissuecols, labels=shortlab, name='Tissue', limits=tissues) +
  scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +
  scale_shape_manual(values=c('1'=21,'0'=24)) +
  
  # fill
  geom_label_repel(data=finemap[!is.na(label)],
                   aes(x=finemap_clpp_global,y=finemap_clpp_local,label=shortl,fill=tissue,colour=tissue),
                   seed=1,
                   alpha=0.4,
                   size=3,
                   show.legend=F) + 
  
  # no fill
  geom_label_repel(data=finemap[!is.na(label)],
                   aes(x=finemap_clpp_global,y=finemap_clpp_local,label=shortl,colour=tissue),
                   seed=1,
                   fill=NA,
                   size=3,
                   show.legend=F) + 
  
  # black text 
  geom_label_repel(data=finemap[!is.na(label)],
                   aes(x=finemap_clpp_global,y=finemap_clpp_local,label=shortl),
                   seed=1,
                   fill=NA,
                   colour='black',
                   label.size=NA,
                   segment.size=NA,
                   size=3,
                   show.legend=F) + 
  
  labs(x='GlobalAA FINEMAP CLPP', y='LocalAA FINEMAP CLPP') +
  theme_classic() +
  #coord_trans(x='log10',y='log10',xlim=c(0.0008,1), ylim=c(0.0008,1), expand=c(0,0), ) +
  coord_cartesian(xlim=c(0.0008, 1), ylim=c(0.0008, 1)) +
  scale_y_continuous(expand = c(0,0), trans='log10', breaks = c(.001, .01, .1, 1), labels=c('0.001','0.01','0.1','1')) +
  scale_x_continuous(expand = c(0,0), trans='log10', breaks = c(.001, .01, .1, 1), labels=c('0.001','0.01','0.1','1')) +
  guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=3),ncol=2),fill=FALSE,label=F,
         shape = F) + 
  theme(panel.grid=element_blank(),
        legend.position='none'
        #legend.position=c(0.8,0.8),
        #legend.key.height = unit(1.5, 'lines'),
        #legend.justification='center',
        #legend.title=element_blank(),
        #legend.margin=margin(b = -2, unit='mm'),
        #legend.direction='horizontal'
  ) +
  geom_abline(linetype='dashed') 

png(paste0(plot_dir,'/Fig4-FINEMAP.png'), width=5,height=5,unit='in',res=300)
g2
dev.off()

# plot together
pdf(paste0(plot_dir,'/Fig4-both-coloc.pdf'), width=8,height=4)
grid.arrange(g, g2, layout_matrix=rbind(c(1,1,2,2),c(1,1,2,2)))
dev.off()


# make legend for labels 
sub = master_coloc[!is.na(shortl)]
sub = sub[order(shortl)]
sub[,y := rev(1:6)]
sub[,x := 1:6]
pdf(paste0(plot_dir,'/Fig4-text-legend.pdf'), width=10,height=4)
ggplot(sub, aes(x=x, y=y, colour=tissue)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim=c(6,30)) +
  geom_label_repel(aes(label=shortl, colour=tissue, fill=tissue),
                   nudge_x=8-sub[,x],
                   segment.size = NA,
                   seed=1, 
                   alpha=0.4,
                   show.legend = F) +
  geom_label_repel(aes(label=shortl, colour=tissue),
                   nudge_x=8-sub[,x],
                   segment.size = NA,
                   seed=1, 
                   alpha=0.4,
                   show.legend = F,
                   fill=NA) +
  geom_label_repel(aes(label=shortl),
                   nudge_x=8-sub[,x],
                   segment.size = NA,
                   seed=1, 
                   label.size=NA,
                   colour='black',
                   show.legend = F,
                   fill=NA) +
  scale_colour_manual(values=tissuecols, labels=shortlab, name='Tissue', guide='none') +
  scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +
  geom_text(aes(label=label,x=10), colour='black', hjust=0)
dev.off()

  
  

#################################################################################################

# fst regression
# TRANSFERRED 
load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/lmer_coef-to-plot.RData')
ann_text = data.frame(coefficient = c('methodlocal:max_within_fst','max_within_fst'),
                      method = c('COLOC PP4', 'COLOC PP4'),
                      estimate = c(0.38,-0.6),
                      label = c('**','***'))

g = ggplot(lmer_coef, aes(y = coefficient, x = estimate)) +
  geom_vline(xintercept=0, linetype='dashed') +
  geom_point() +
  geom_errorbar(aes(xmin=estimate-se, xmax=estimate+se), width=.2) +
  theme_classic() +
  facet_wrap(~method, ncol=2) + 
  theme(axis.text.y=element_text(colour='black'),
        axis.title.y=element_blank()) +
  labs(x='Coefficient estimate') +
  scale_y_discrete(limits = rev(c('methodlocal',
                                  'AFR_EUR_fst',
                                  'max_within_fst',
                                  'methodlocal:AFR_EUR_fst',
                                  'methodlocal:max_within_fst')),
                   labels=c('methodlocal:max_within_fst'=expression(LocalAA:F[ST*","*within]),
                            'methodlocal:AFR_EUR_fst'=expression(LocalAA:F[ST*","*between]),
                            'methodlocal'='LocalAA',
                            'max_within_fst'=expression(F[ST*","*within]),
                            'AFR_EUR_fst'=expression(F[ST*","*between]))) +
  geom_text(data = ann_text, aes(label = label), vjust=0.8, colour='red')
