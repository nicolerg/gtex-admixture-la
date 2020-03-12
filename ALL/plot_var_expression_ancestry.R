# plot variance explained in gene expression by global and local ancestry (see variance_explained_ancestry_gene_expression)

library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

methodcols <- c(global='#81c0c7',LAVA='#CC7EAB')

tissuelabs <- c(Adipose_Subcutaneous="Subc. adipose",
	Artery_Tibial="Tibial artery",
	Lung="Lung",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Skin_Not_Sun_Exposed_Suprapubic="NSE skin")

tissuecols <- c(Adipose_Subcutaneous="#FF6600",
	Adipose_Visceral_Omentum="#FFAA00",
	Adrenal_Gland="#33DD33",
	Artery_Aorta="#FF5555",
	Artery_Coronary="#FFAA99",
	Artery_Tibial="#FF0000",
	Brain_Amygdala="#EEEE00",
	Brain_Anterior_cingulate_cortex_BA24="#EEEE00",
	Brain_Caudate_basal_ganglia="#EEEE00",
	Brain_Cerebellar_Hemisphere="#EEEE00",
	Brain_Cerebellum="#EEEE00",
	Brain_Cortex="#EEEE00",
	Brain_Frontal_Cortex_BA9="#EEEE00",
	Brain_Hippocampus="#EEEE00",
	Brain_Hypothalamus="#EEEE00",
	Brain_Nucleus_accumbens_basal_ganglia="#EEEE00",
	Brain_Putamen_basal_ganglia="#EEEE00",
	"Brain_Spinal_cord_cervical_c-1"="#EEEE00",
	Brain_Substantia_nigra="#EEEE00",
	Breast_Mammary_Tissue="#33CCCC",
	"Cells_EBV-transformed_lymphocytes"="#CC66FF",
	Cells_Cultured_fibroblasts="#AAEEFF",
	Colon_Sigmoid="#EEBB77",
	Colon_Transverse="#CC9955",
	Esophagus_Gastroesophageal_Junction="#8B7355",
	Esophagus_Mucosa="#552200",
	Esophagus_Muscularis="#BB9988",
	Heart_Atrial_Appendage="#9900FF",
	Heart_Left_Ventricle="#660099",
	Liver="#AABB66",
	Lung="#99FF00",
	Kidney_Cortex="yellow",
	Minor_Salivary_Gland="#99BB88",
	Muscle_Skeletal="#AAAAFF",
	Nerve_Tibial="#FFD700",
	Ovary="#FFAAFF",
	Pancreas="#995522",
	Pituitary="#AAFF99",
	Prostate="#DDDDDD",
	Skin_Not_Sun_Exposed_Suprapubic="#0000FF",
	Skin_Sun_Exposed_Lower_leg="#7777FF",
	Small_Intestine_Terminal_Ileum="#555522",
	Spleen="#778855",
	Stomach="#FFDD99",
	Testis="#AAAAAA",
	Thyroid="#006600",
	Uterus="#FF66FF",
	Vagina="#FF5599",
	Whole_Blood="#FF00BB")

i <-1 
dt_list <- list()
for (t in tissues){
	dt <- fread(sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/%s_var_explained_by_ancestry.tsv',t),sep='\t',header=T)
	dt[,tissue := t]
	dt_list[[i]] <- dt
	i <- i + 1
}

merged <- rbindlist(dt_list)
head(merged)

merged_melt <- melt(merged, id.vars=c('tissue','gene'))
head(merged_melt)

pdf('~/gtex-admix/plots/august/var-explained-expression-ancestry.pdf',width=12,height=6)
ggplot(merged_melt, aes(value)) +
	geom_density(aes(colour=factor(variable),fill=factor(variable)), alpha=0.5) +
	facet_wrap(~tissue,ncol=3) + 
	theme_bw() + 
	labs(x='Variance explained in gene expression',y='Density') +
	scale_colour_manual(values=c(variance_explained_la=methodcols[['LAVA']],variance_explained_ga=methodcols[['global']]),guide='none') +
	scale_fill_manual(values=c(variance_explained_la=methodcols[['LAVA']],variance_explained_ga=methodcols[['global']]), labels=c(variance_explained_la='Local ancestry',variance_explained_ga='Global ancestry')) +
	theme(legend.title=element_blank())
dev.off()

pdf('~/gtex-admix/plots/august/var-explained-expression-ancestry-scatter.pdf',width=12,height=6)
ggplot(merged, aes(x=variance_explained_ga, y=variance_explained_la)) +
	geom_abline(linetype='dashed') +
	geom_point(alpha=0.5) +
	facet_wrap(~tissue,ncol=3) + 
	theme_bw() + 
	labs(x='Var. explained in gene expression by GA',y='Var. explained in gene expression by LA') 
dev.off()

# shuffle rows
merged <- merged[sample(1:nrow(merged),nrow(merged),replace=F)]

pdf('~/gtex-admix/plots/august/var-explained-expression-ancestry-scatter-combined.pdf',width=8,height=6)
ggplot(merged, aes(x=variance_explained_ga, y=variance_explained_la, colour=factor(tissue),label=gene)) +
	geom_abline(linetype='dashed') +
	geom_point(alpha=0.5) +
	theme_bw() + 
	labs(x='Var. explained in gene expression by GA',y='Var. explained in gene expression by LA') +
	scale_colour_manual(values=tissuecols, labels=tissuelabs) +
	theme(legend.title = element_blank()) +
	geom_text_repel(data=merged[variance_explained_la > 0.157],
		nudge_x = 0.1 - merged[variance_explained_la > 0.157, variance_explained_ga],
		direction = 'y',
		hjust=0)
dev.off()

# look at correlation between each PC and each local ancestry variable 

i <-1 
dt_list <- list()
for (t in tissues){
	dt <- fread(sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/%s_var_ancestry_correlation.tsv',t),sep='\t',header=T)
	dt[,tissue := t]
	dt_list[[i]] <- dt
	i <- i + 1
}

corr <- rbindlist(dt_list)
head(corr)

corr <- unique(corr, by='gene')
corr_melted <- melt(corr, id.var=c('gene','tissue'))

print(head(corr_melted))

pdf('~/gtex-admix/plots/august/ancestry-cor-per-gene.pdf',width=8,height=6)
ggplot(corr_melted, aes(value)) +
	geom_histogram(bins=50) +
	facet_wrap(~variable) +
	geom_vline(xintercept=c(-1,1,0), linetype='dashed') +
	theme_bw()
dev.off()

corr_melted[,x := sapply(variable, function(j) unname(unlist(strsplit(as.character(j), '_')))[1])]
corr_melted[,x := paste0('g',x)]
corr_melted[,group := sapply(variable, function(j) unname(unlist(strsplit(as.character(j), '_')))[2])]

pdf('~/gtex-admix/plots/august/ancestry-cor-per-gene.pdf',width=5,height=3)
ggplot(corr_melted, aes(x=factor(x), y=value, fill=factor(group))) + 
	geom_hline(yintercept=0,linetype='dashed') +
	geom_boxplot(colour='black') + 
	theme_bw() + 
	scale_fill_manual(values=c(ASN="#ADFF33", AFR="#FF9900")) +
	labs(x='Genotype PC',y='Correlation between LA at TSS and gPC') +
	theme(axis.title.x=element_blank(),
		axis.text.x=element_text(colour='black'),
		legend.title=element_blank(),
		legend.direction='horizontal',
		legend.position=c(0.85,0.9),
		panel.border = element_rect(colour = "black", fill=NA, size=1),
		legend.box.background = element_rect(colour = "black", size=1)) +
	ylim(c(-1,1))
dev.off()
