#!/bin/R

library(data.table)
library(ggplot2)
library(ggcorrplot)

# ADMIXTURE global ancestry estimates
admix <- fread('merged_hg38_filt_MAF10_nomissing_pruned.3.Q.subj',header=F)
colnames(admix) <- c('SUBJID','SUB_POP','SUPER_POP','AFR','EUR','ASN')
# GTEx genotype PCs
pcs <- fread('GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',sep='\t',header=T)
# global ancestry estimates from RFMix
rfmix <- fread('gtex_ai_global_all.txt',sep='\t',header=T)
# GTEx subject phenotypes 
gtex_meta <- fread('GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
gtex_meta <- gtex_meta[,.(SUBJID, RACE)]

# merge RFMix and ADMIXTURE
rfmix_melt = melt(rfmix, id.vars=c('ID'), variable.name = 'GA', value.name = 'RFMix')
admix_melt = melt(admix, id.vars=c('SUBJID','SUB_POP','SUPER_POP'), variable.name = 'GA', value.name = 'ADMIXTURE')
m_melt = merge(admix_melt, rfmix_melt, by.x=c('SUBJID','GA'),by.y=c('ID','GA'))

# plot RFMix vs ADMIXTURE global ancestry estimates for each superpopulation
pdf('~/Desktop/admixture-v-rfmix.pdf',width=8,height=3)
ggplot(m_melt, aes(x=RFMix, y=ADMIXTURE)) +
  geom_abline(linetype='dashed') +
  geom_point() +
  scale_color_identity() +
  theme_classic() +
  facet_wrap(~GA) +
  xlim(c(0,1)) +
  ylim(c(0,1))
dev.off()

# merge gPCs with subject metadata
pcs[,SUBJID:=sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
pcs[,c('FID','IID'):=NULL]
m <- merge(gtex_meta, pcs, by='SUBJID')
# merge with RFMix
colnames(rfmix) <- paste0(colnames(rfmix),'_rfmix')
m <- merge(m, rfmix, by.x='SUBJID',by.y='ID_rfmix')
# merge with ADMIXTURE
admix_gtex <- admix[grepl('GTEX',SUBJID)]
colnames(admix_gtex) = c('SUBJID','SUB_POP','SUPER_POP','AFR_admixture','EUR_admixture','ASN_admixture')
m <- merge(m, admix_gtex, by='SUBJID')

# plot ADMIXTURE global ancestry proportions; label individuals with GTEx-reported ancestry 
m <- m[order(EUR_admixture)]
melt <- melt(m, measure.vars=colnames(m)[grepl('admixture',colnames(m))])
melt[,variable := factor(variable, levels=c('EUR_admixture','AFR_admixture','ASN_admixture','EAS_admixture','SAS_admixture'))]
pdf('~/Desktop/admixture-3-superpop-mixed.pdf',width=8,height=3.5)
ggplot(melt, aes(x=SUBJID,fill=variable)) +
  geom_bar(stat='identity',aes(y=value)) +
  theme_classic() +
  scale_fill_manual(values=c(EUR_admixture='blue',
                             AFR_admixture='orange',
                             ASN_admixture='green3',
                             EAS_admixture='green3',
                             SAS_admixture='purple'),
                    name='ADMIXTURE population',
                    labels=c(EUR_admixture='EUR',
                             AFR_admixture='AFR',
                             ASN_admixture='ASN',
                             EAS_admixture='EAS',
                             SAS_admixture='SAS')) +
  coord_cartesian(expand=FALSE, ylim=c(-0.04,1.01)) +
  scale_x_discrete(limits=m[,SUBJID]) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        line = element_blank(),
        legend.position = 'bottom',
        legend.margin = margin(t = -0.4, b=-0.1, unit='cm'),
        legend.spacing = unit(2, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.text = element_text(margin = margin(r = 0.2, l=-0.1, unit = 'cm')),
        legend.key.size = unit(0.5, "lines"),
        axis.title.x = element_blank()) +
  labs(y='ADMIXTURE-reported global ancestry') +
  geom_point(aes(y=-0.02,colour=factor(RACE)),size=0.8) +
  scale_colour_manual(values=c('1'="green3", '2'="orange", '3'="blue", '4'='#FF3D3D', '99'='gray'),
                      labels=c('1'="ASN", '2'="AFR", '3'="EUR", '4'='Native Am.', '99'='Unknown'),
                      breaks=c('3','2','1','99','4'),name='GTEx self-reported ancestry') +
  guides(fill=guide_legend(override.aes = list(colour=NA, shape=15),
                           order=1,
                           title.position="top", title.hjust = 0.5),
         colour=guide_legend(override.aes = list(size=2),
                             order=2,
                             title.position="top", title.hjust = 0.5))
dev.off()

# plot correlations between genotype PCs, ADMIXTURE estimates, and RFMix estimates
m[,SUB_POP := NULL]
m[,SUPER_POP := NULL]
m[,UNK_rfmix := NULL]
colnames(m)[colnames(m)=='PC1'] = 'gPC1'
colnames(m)[colnames(m)=='PC2'] = 'gPC2'
pdf('~/Desktop/admixture-3-superpop-mixed-corr.pdf',width=6,height=6)
ggcorrplot(cor(m[,c(3:4,23:ncol(m)),with=F]),lab=T,hc.order=T)
dev.off()
