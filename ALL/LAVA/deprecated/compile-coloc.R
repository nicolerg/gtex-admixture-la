library(data.table)
library(ggplot2)
library(ggrepel)

# DEPRECATED

cutoff <- 1e-6


# # # load uniq eGenes
# # load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/uniq_egene_',cutoff,'.RData'))
# # # master_uniq_egene



# # for diff lead ============================

# # # load same eGenes, different lead SNP
# # load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_',cutoff,'.RData'))
# # # master_diff_lead

# # # get one SNP per gene
# # master_diff_lead <- master_diff_lead[,.(gene_id, variant_id, pval_nominal, tissue, method)]
# # master_diff_lead <- unique(master_diff_lead, by=c('gene_id','tissue','method'))
# # head(master_diff_lead)
# # g <- master_diff_lead[method=='global',.(gene_id, variant_id, pval_nominal, tissue)]
# # l <- master_diff_lead[method=='LAVA',.(gene_id, variant_id, pval_nominal, tissue)]
# # diff_lead <- merge(g,l,by=c('gene_id','tissue'),suffixes=c('_global','_local'))
# # diff_lead[,pval_diff := -log10(pval_nominal_local) - -log10(pval_nominal_global)]
# # nrow(diff_lead)

# # # merge with coloc results
# # f <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master-coloc-v2.tsv.gz'
# # coloc <- fread(cmd=sprintf("zcat %s", f), sep='\t', header=TRUE)

# # diff_lead <- merge(diff_lead, coloc, by=c('gene_id', 'tissue'))
# # diff_lead[,p4 := as.numeric(p4)]
# # diff_lead[,p3 := as.numeric(p3)]
# # nrow(diff_lead)
# # save(diff_lead, file=paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/diff_lead-',cutoff,'.RData'))
# load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/diff_lead-',cutoff,'.RData'))

# master_coloc <- diff_lead[order(pval_diff, decreasing=TRUE)]

# # 1 p4 value per gene per trait per tissue

# # filter down to p4 > 0.01
# nrow(master_coloc)
# # master_coloc <- master_coloc[p4 > 0.01]
# # nrow(master_coloc)

# # p4
# master_coloc[,colour := ifelse(p4 > 0.5, 'red','black')]
# master_coloc[,label := ifelse(p4 > 0.5, trait, NA)]
# jpeg('~/gtex-admix/plots/final_figures/v8-p4-coloc.jpg',res=150,width=1800,height=1000)
# ggplot(master_coloc, aes(x=pval_diff, y=p4, colour=colour, label=label)) +
# 	geom_point() +
# 	geom_text_repel(size=3) +
# 	facet_wrap(~tissue, ncol=3) +
# 	theme_bw() +
# 	scale_colour_manual(values=c(black='black',red='red'),guide='none') +
# 	geom_vline(xintercept=0,linetype='dashed') +
# 	labs(x='-log10(Local pvalue) - -log10(Global pvalue)', y='Colocalization p-value (p4)')
# dev.off()

# # p3 - probably more informative
# master_coloc[,colour := ifelse(p3 > 0.8, 'red','black')]
# jpeg('~/gtex-admix/plots/final_figures/v8-p3-coloc.jpg',res=150,width=1800,height=1000)
# ggplot(master_coloc, aes(x=pval_diff, y=p3, colour=colour)) +
# 	geom_point() +
# 	facet_wrap(~tissue, ncol=3) +
# 	theme_bw() +
# 	scale_colour_manual(values=c(black='black',red='red'),guide='none') +
# 	geom_vline(xintercept=0,linetype='dashed',colour='blue') +
# 	labs(x='-log10(Local pvalue) - -log10(Global pvalue)', y='Colocalization p-value (p3)')
# dev.off()

# # filter down to p3 > 0.8
# coloc_filt <- master_coloc[p3 > 0.8]
# nrow(coloc_filt)

# jpeg('~/gtex-admix/plots/final_figures/v8-p3-filt-coloc.jpg',res=150,width=1800,height=1000)
# ggplot(coloc_filt, aes(x=pval_diff, y=p3, colour=gene_id)) +
# 	geom_point() +
# 	facet_wrap(~tissue, ncol=3) +
# 	theme_bw() +
# 	geom_vline(xintercept=0,linetype='dashed') +
# 	labs(x='-log10(Local pvalue) - -log10(Global pvalue)', y='Colocalization p-value (p3)') +
# 	theme(legend.position = 'none') 
# dev.off()

# counts <- data.table(table(coloc_filt[,gene_id],coloc_filt[,tissue]))
# coloc_filt <- coloc_filt[,.(gene_id, tissue, pval_diff)] # all p3 > 0.8
# coloc_filt <- unique(coloc_filt)
# coloc_filt <- merge(coloc_filt, counts, by.x=c('gene_id','tissue'), by.y=c('V1','V2'))
# head(coloc_filt)

# #write.table(unique(coloc_filt[,gene_id]), '~/gtex-admix/metadata/lead-diff-gene-id.txt', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
# gene_labels <- fread('~/gtex-admix/metadata/lead-diff-gene-id-name-map.txt',sep='\t',header=TRUE)
# head(gene_labels)
# coloc_filt <- merge(coloc_filt, gene_labels, by.x='gene_id', by.y='Gene stable ID version')

# coloc_filt[,label := ifelse(N > 18, `Gene name`, NA)]
# jpeg('~/gtex-admix/plots/final_figures/v8-high-p3-coloc-bygene.jpg',res=150,width=1800,height=1000)
# ggplot(coloc_filt, aes(x=pval_diff, y=N, colour=gene_id, label=label)) +
# 	geom_point() +
# 	geom_text_repel(size=3, xlim=c(2,NA)) + 
# 	facet_wrap(~tissue, ncol=3) +
# 	theme_bw() +
# 	geom_vline(xintercept=0,linetype='dashed') +
# 	labs(x='-log10(Local pvalue) - -log10(Global pvalue)', y='Number of GWAS traits where eGene has p3 > 0.8') +
# 	theme(legend.position = 'none') 
# dev.off()

# look at when there is a big difference in pvalue - when Global is more significant, are they more polygenic traits?




##############################################################################
## Look at my own coloc results 

library(data.table)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

tissuecols <- c(Adipose_Subcutaneous="#FF6600",
	Artery_Tibial="#FF0000",
	Lung="#99FF00",
	Muscle_Skeletal="#AAAAFF",
	Nerve_Tibial="#FFD700",
	Skin_Not_Sun_Exposed_Suprapubic="#0000FF")

tissuelabs <- c(Adipose_Subcutaneous="Adipose",
	Artery_Tibial="Tibial artery",
	Lung="Lung",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Skin_Not_Sun_Exposed_Suprapubic="Skin (NSE)")

# # get subset of genes we're interested in 

# load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/diff_lead-',cutoff,'.RData')) # same eGene, different lead SNP 
# dtlist <- list()
# i <- 1
# for (tissue in tissues){
# 	#gene_list <- data.table(gene_id=unique(diff_lead[tissue==tissue & (p3 > 0.8 | p4 > 0.5),gene_id]))
# 	gene_list <- diff_lead[tissue==tissue & (p3 > 0.8 | p4 > 0.5)]
# 	colnames(gene_list)[colnames(gene_list)=='p3'] <- 'gtex_v8_p3'
# 	colnames(gene_list)[colnames(gene_list)=='p4'] <- 'gtex_v8_p4'
# 	# trait_counts <- data.table(table(diff_lead[tissue==tissue & (p3 > 0.8 | p4 > 0.5),trait]))
# 	# write.table(trait_counts, paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/genelist/',tissue,'.traits.txt'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

# 	# # filter allpairs 

# 	# write.table(gene_list, paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/genelist/',tissue,'.genelist.txt'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
# 	# system(sprintf('python ~/gtex-admix/scripts/LAVA/annotation/parse_allpairs_for_coloc_smr.py %s &',tissue))
# 	dtlist[[i]] <- gene_list
# 	i <- i+1
# }

# gene_subset <- data.table(rbindlist(dtlist))

# read in MY coloc results

get_ref_snp <- function(x){
	return(paste(unname(unlist(strsplit(x,'_')))[1:2],collapse='_'))
}

# load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_',cutoff,'.RData'))
# head(master_diff_lead)
# print(length(unique(master_diff_lead[,gene_id])))
# master_diff_lead[method=='global',ref_snp_global := sapply(variant_id, get_ref_snp)]
# master_diff_lead[method=='LAVA',ref_snp_local := sapply(variant_id, get_ref_snp)]
# master_diff_lead[,variant_id := NULL]
# gene_subset <- master_diff_lead

my_results <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/output/gtex-admixed-difflead/2019-05-02_15-00-25.673667_117admixed_difflead/merged_coloc.txt'

#coloc <- fread(cmd=sprintf('zcat %s',my_results), header=TRUE, sep='\t')
coloc <- fread(my_results, header=TRUE, sep='\t')

# format COLOC results to be mergable with gene_subset

coloc[,ref_snp := paste0('chr',ref_snp)]

# gene_subset[,ref_snp_global := sapply(variant_id_global, get_ref_snp)]
# gene_subset[,ref_snp_local := sapply(variant_id_local, get_ref_snp)]
# gene_subset[,variant_id_local := NULL]
# gene_subset[,variant_id_global := NULL]

# get tissue from COLOC Adipose_Subcutaneous_filtered_af_global_allpairs_mike_tsv_gz
coloc[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]

# get method
coloc[, method := ifelse(grepl('global', eqtl_file),'global','LAVA')] 
coloc[, eqtl_file := NULL]

coloc[, gwas_trait := gsub('coloc_imputed_','',gwas_trait)]
coloc[, gwas_trait := gsub('\\..*','',gwas_trait)]
coloc[, base_gwas_file := NULL]

# # look at correlation between my global results and the GTEx v8 results

# global_only <- merge(coloc[method=='global'], gene_subset, by.x=c('gwas_trait','feature','tissue'), by.y=c('trait','gene_id','tissue'))
# global_only <- unique(global_only)
# nrow(global_only)
# head(global_only)

# pdf('~/gtex-admix/plots/final_figures/compare_coloc_p-values_difflead.pdf',height=4,width=4)
# ggplot(global_only) +
# 	geom_point(aes(x=gtex_v8_p3, y=clpp_h3),alpha=0.5) +
# 	labs(x='GWAS subgroup p3', y='Admixed p3') +
# 	geom_abline(colour='red', linetype='dashed') +
# 	theme_bw()

# ggplot(global_only) +
# 	geom_point(aes(x=gtex_v8_p4, y=clpp_h4),alpha=0.5) +
# 	labs(x='GWAS subgroup p4', y='Admixed p4') +
# 	geom_abline(colour='red', linetype='dashed') +
# 	theme_bw()
# dev.off()

# now compare only MY coloc results between global and local 

coloc[, n_snps:=NULL]

# ref_snp from coloc results is the lead SNP for that eGene

# genes_local <- gene_subset[,.(gene_id, tissue, trait, ref_snp_local)]
# genes_global <- gene_subset[,.(gene_id, tissue, trait, ref_snp_global)]

# load(paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/diff_lead-',cutoff,'.RData'))

coloc_global <- unique(coloc[method=='global'])
coloc_global[,method:=NULL]
# coloc_global <- unique(merge(coloc_global, genes_global, by.x=c('feature','tissue','gwas_trait'), by.y=c('gene_id','tissue','trait'),all.y=TRUE))
# colnames(coloc_global)[colnames(coloc_global)=='ref_snp_global'] <- 'ref_snp_tied'

coloc_local <- unique(coloc[method=='LAVA'])
coloc_local[,method:=NULL]
# coloc_local <- unique(merge(coloc_local, genes_local, by.x=c('feature','tissue','gwas_trait'), by.y=c('gene_id','tissue','trait'),all.y=TRUE))
# colnames(coloc_local)[colnames(coloc_local)=='ref_snp_local'] <- 'ref_snp_tied'

nrow(coloc_global)
nrow(coloc_local)

m <- merge(coloc_local, coloc_global, by=c('feature','tissue','gwas_trait','ref_snp'), suffixes=c('_local','_global'))
head(m)
nrow(m)

# filter by gene list to remove extra genes from Skin (?)

difflead_genes <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/all_difflead/difflead_gene_list.txt',sep='\t',header=FALSE)
head(difflead_genes)

m <- m[feature %in% difflead_genes[,V1]]

# # add p-value difference column
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-06.RData')

# diff_lead <- master_diff_lead[,.(gene_id, variant_id, pval_nominal, method, tissue)]
# diff_lead[,ref_snp := sapply(variant_id, get_ref_snp)]

# m <- merge(coloc, diff_lead, by.x=c('feature','tissue','method'),by.y=c('gene_id','tissue','method'),suffixes=c('_gwas','_eqtl'))
# m <- unique(m, by=c('feature','tissue','method','gwas_trait'))

# g <- m[method=='global']
# g[,method:=NULL]
# l <- m[method=='LAVA']
# l[,method:=NULL]

# m <- merge(g, l, by=c('feature','tissue','gwas_trait'),suffixes=c('_global','_local'))

# nrow(m)
# master_coloc <- m[complete.cases(m)]
# nrow(master_coloc)
# head(master_coloc)

# look to see if local lead SNP has higher coloc prob than global lead SNP when local lead SNP is more significant

master_coloc <- m
# master_coloc[,pval_diff := -log10(pval_nominal_local) + log10(pval_nominal_global)]
# master_coloc[,line_colour := ifelse(clpp_h4_global > clpp_h4_local, 'global', 'LAVA')]

methodcols <- c(global='#8DD3C7',LAVA='#FDB462')
methodlabs <- c(global='Global',LAVA='Local')

# pdf('~/gtex-admix/plots/final_figures/COLOC-local-global-p4-change-difflead.pdf', width=10,height=4)
# ggplot(master_coloc) + 
# 	geom_segment(data=master_coloc[abs(clpp_h4_global - clpp_h4_local) < 0.05], aes(x=pval_diff, xend=pval_diff, y=clpp_h4_global, yend=clpp_h4_local), colour='gray') +
# 	geom_point(data=master_coloc[abs(clpp_h4_global - clpp_h4_local) < 0.05], aes(x=pval_diff, y=clpp_h4_global), color='gray') +
# 	geom_point(data=master_coloc[abs(clpp_h4_global - clpp_h4_local) < 0.05], aes(x=pval_diff, y=clpp_h4_local), color='gray') +
# 	geom_segment(data=master_coloc[abs(clpp_h4_global - clpp_h4_local) > 0.05], aes(x=pval_diff, xend=pval_diff, y=clpp_h4_global, yend=clpp_h4_local, colour=line_colour)) +
# 	geom_point(data=master_coloc[abs(clpp_h4_global - clpp_h4_local) > 0.05], aes(x=pval_diff, y=clpp_h4_global), color='#8DD3C7') +
# 	geom_point(data=master_coloc[abs(clpp_h4_global - clpp_h4_local) > 0.05], aes(x=pval_diff, y=clpp_h4_local), color='#FDB462') +
# 	scale_colour_manual(values=methodcols, labels=methodlabs, name='Ancestry') +
# 	theme_bw() +
# 	labs(x='-log10(Local p-value) - -log10(Global p-value)',y='Colocalization prob. (p4)')
# dev.off()

master_coloc[, label:= ifelse((clpp_h4_global < 0.70 & clpp_h4_local > 0.8) | (clpp_h4_local < 0.70 & clpp_h4_global > 0.8), paste0(feature,':',gwas_trait), NA)]
# master_coloc[clpp_h4_global < 0.3 & clpp_h4_local > 0.6]

# master_coloc[, label := NA]
# master_coloc[feature == 'ENSG00000229391.7', label := paste0('HLA-DRB6:',gwas_trait)]
# master_coloc[feature == 'ENSG00000130649.9', label := paste0('CYP2E1:',gwas_trait)]


# pdf('~/gtex-admix/plots/final_figures/COLOC-local-global-p4-scatter-difflead.pdf', width=6,height=4)
# ggplot(master_coloc, aes(x=clpp_h4_global, y=clpp_h4_local)) +
# 	geom_point(aes(x=clpp_h4_global, y=clpp_h4_local, colour=pval_diff),alpha=0.7) +
# 	scale_colour_gradient2(low='#8DD3C7', high='#FDB462', mid='gray', name='-log10(Local p-val) -\n-log10(Global p-val)') +
# 	labs(x='Global p4', y='Local p4') +
# 	theme_bw() +
# 	geom_abline(colour='red',linetype='dashed') +
# 	geom_text_repel(data=master_coloc[clpp_h4_global < 0.3 & clpp_h4_local > 0.6],aes(label=label),size=3.5,ylim=c(1,1.6),direction='y',force=100) +
# 	scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
# 	scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),limits=c(0,1.6))

# dev.off()

jpeg('~/gtex-admix/plots/final_figures/COLOC-local-global-p4-scatter-bytissue.jpg', width=1200,height=800,res=150)
ggplot(master_coloc, aes(x=clpp_h4_global, y=clpp_h4_local)) +
	geom_point(aes(x=clpp_h4_global, y=clpp_h4_local, colour=tissue),alpha=0.7) +
	#scale_colour_gradient2(low='#8DD3C7', high='#FDB462', mid='gray', name='-log10(Local p-val) -\n-log10(Global p-val)') +
	scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue') +
	labs(x='Global p4', y='Local p4') +
	theme_bw() +
	geom_abline(colour='red',linetype='dashed') 
dev.off()

jpeg('~/gtex-admix/plots/final_figures/COLOC-local-global-p4-scatter-bytissue-labelled.jpg', width=1600,height=800,res=150)
ggplot(master_coloc, aes(x=clpp_h4_global, y=clpp_h4_local)) +
	geom_point(aes(x=clpp_h4_global, y=clpp_h4_local, colour=tissue),alpha=0.7) +
	scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue') +
	xlim(c(0,3)) +
	ylim(c(0,1.5)) +
	geom_text_repel(aes(label=label, colour=tissue),size=3, xlim=c(1.1,NA), direction='y', hjust=0) +
	labs(x='Global p4', y='Local p4') +
	theme_bw() +
	geom_abline(colour='red',linetype='dashed') 

dev.off()


# TRY TO AVOID BELOW CODE IF POSSIBLE =======================================================

# filt_p4 <- filt[max_p4 > 0.8]
# head(filt_p4)
# filt_p4[diff < -1.5]
# q()

# t.test(filt[diff > 0, max_p4], filt[diff < 0, max_p4])
# cutoffs <- seq(0,3,0.1)
# pvals <- rep(1, length(cutoffs))
# N <- rep(1, length(cutoffs))
# i <- 1
# for (cutoff in seq(0,3,0.1)){
# 	N[i] <- length(filt[diff > cutoff, max_p4]) + length(filt[diff < -cutoff, max_p4])
# 	pvals[i] <- t.test(filt[diff > cutoff, max_p4], filt[diff < -cutoff, max_p4])$p.value
# 	i <- i+1
# }

# ttest <- data.table(diff_cutoff=cutoffs, pvalue=pvals, N=N)

# # filt[!(diff > 1.5 & max_p4 > 106), trait := NA]
# # filt[!is.na(trait)]

# # jpeg("~/gtex-admix/plots/combined-coloc-scatter.jpg", width=1800, height=600, res=200)
# # ggplot() +
# # 	geom_point(data=filt, aes(x=diff, y=max_p4)) +
# # 	theme_bw() +
# # 	theme(axis.text.x=element_blank(),
# # 		axis.title.x=element_blank()) + 
# # 	labs(y='Max p4 p-value (-log10)') +
# # 	geom_text(data=filt, aes(x=diff+0.1, y=max_p4, label=trait), size=3, hjust=0) +
# # 	geom_vline(xintercept=0, linetype='dashed') +
# # 	xlim(c(-5,5.5))
# # dev.off()

# jpeg("~/gtex-admix/plots/combined-coloc-scatter-unlabelled.jpg", width=1800, height=600, res=200)
# ggplot() +
# 	geom_point(data=filt, aes(x=diff, y=max_p4)) +
# 	theme_bw() +
# 	theme(axis.title.x=element_blank()) + 
# 	labs(y='Max p4 value') +
# 	geom_vline(xintercept=0, linetype='dashed') 
# dev.off()

# jpeg("~/gtex-admix/plots/combined-coloc-scatter-filt.jpg", width=1800, height=600, res=200)
# ggplot() +
# 	geom_point(data=filt_p4, aes(x=diff, y=max_p4)) +
# 	theme_bw() +
# 	theme(axis.title.x=element_blank()) + 
# 	labs(y='Max p4 value') +
# 	geom_vline(xintercept=0, linetype='dashed') 
# dev.off()

# pdf("~/gtex-admix/plots/combined-coloc.pdf", width=10, height=4)

# ggplot(ttest) +
# 	geom_point(aes(x=diff_cutoff, y=-log10(pvalue)), colour='red') +
# 	geom_text(aes(x=diff_cutoff, y=-log10(pvalue)+0.05, label=N), angle=90, hjust=0) +
# 	labs(x='P-value difference filter',y='t-test p-value (-log10)') +
# 	theme_bw() +
# 	geom_hline(yintercept=-log10(0.05), linetype='dashed')

# ggplot() +
# 	geom_histogram(data=filt[1:500], aes(x=max_p4, fill='pink'),alpha=0.5,bins=100) +
# 	geom_histogram(data=tail(filt,500), aes(x=max_p4, fill='lightblue'),alpha=0.5,bins=100) +
# 	scale_fill_manual(values=c(pink='red',lightblue='blue'),
# 		labels=c(
# 			pink='Top 500 diff in LAVA',
# 			lightblue='Top 500 diff in Global'),
# 		name='Max p4 per gene and tissue pair') +
# 	labs(x='-log10(pvalue)') +
# 	theme_bw() 

# ggplot() +
# 	geom_histogram(data=filt, aes(x=diff),fill='gray',alpha=0.5,bins=100) +
# 	theme_bw() +
# 	labs(x='-log10(LAVA) - -log10(Global)')

# ggplot() +
# 	geom_histogram(data=filt[1:500], aes(x=diff, fill='red'),alpha=0.5,bins=100) +
# 	geom_histogram(data=tail(filt,500), aes(x=-diff, fill='blue'),alpha=0.5,bins=100) +
# 	scale_fill_manual(values=c(red='red',blue='blue'),
# 		labels=c(
# 			red='Top 500 LAVA',
# 			blue='Top 500 Global'),
# 		name='Difference b/w LAVA and Global -log10(pval)') +
# 	labs(x='Difference in -log10(pvalue)') +
# 	theme_bw() 

# ggplot(data=filt, aes(x=diff, y=max_p4)) +
# 	geom_point() +
# 	labs(x='-log10(LAVA) - -log10(Global)', y='Max p4')

# dev.off()








