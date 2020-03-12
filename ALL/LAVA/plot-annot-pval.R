library(ggplot2)
library(data.table)
require('OneR')
library(reshape)
library(stats)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
tissue <- args[2]
full <- args[3]
allp <- args[4]
caddmap <- args[5]

# infile <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/Lung.diff.SNP.annot.tsv'
# tissue <- 'Lung'
# full <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/Lung.diff.eSNP.annotated.tsv'
# allp <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/Lung-LAVA-global-merged-allpairs-filt.tsv.gz'
# caddmap <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/Lung.cadd.map.txt'

#########################################################################################

# # plot annotation of non-overlapping lead SNPs across bins of differences in p-value

# summ <- fread(infile, sep='\t', header=TRUE)
# summ[,pval_diff := -log10(lava_pval) - -log10(global_pval)]

# # divide data into bins

# bins <- bin(summ[,pval_diff],nbins=10,method='length',labels=1:10)
# summ[,bin := bins]

# bin_labels <- bin(summ[,pval_diff],nbins=10,method='length')
# bin_lab <- data.table(bins=bins, labels=bin_labels)
# bin_lab <- bin_lab[!duplicated(bins)]
# bin_lab <- bin_lab[order(bins)]
# x <- data.table(table(summ[,bin]))
# x <- x[N != 0]
# bin_lab[,counts := x[,N]] # V1, N

# # data <- melt(merged, id=c('gene_id','pval_nominal.global','pval_nominal.lava','pval_diff','bin'))
# # head(data)

# # jpeg(paste0("/users/nicolerg/gtex-admix/plots/",tissue,"-pval-diff-v-annot.jpg"),width=1600,height=1000,res=200)
# # ggplot(data, aes(x=bin,y=value,colour=factor(variable))) + 
# # 	geom_boxplot() +
# # 	labs(x='Difference in nominal p-value of top SNPs (global - LAVA)', y='Observed prob. of annotation in top SNPs') +
# # 	scale_colour_manual(values=c(prob_annot.global='blue',prob_annot.lava='red'), 
# # 		labels=c(prob_annot.global='Global',prob_annot.lava='LAVA'),name='Method') +
# # 	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
# # 	scale_x_discrete(breaks=bin_lab[,bins], labels=bin_lab[,labels]) +
# # 	coord_cartesian(ylim=c(0,1.1)) +
# # 	annotate('text',x=bin_lab[,bins],y=1.1,label=bin_lab[,counts])
# # dev.off()

# # binarize 
# summ[,prob_annot.global := ifelse(global_annot > 0, 1, 0)]
# summ[,prob_annot.lava := ifelse(lava_annot > 0, 1, 0)]

# d <- data.table(bin=unique(summ[,bin]),pct_lava=0,pct_global=0)

# for (b in unique(summ[,bin])){
# 	sub <- summ[bin == b]
# 	print(sum(sub[,prob_annot.lava])/nrow(sub))
# 	d[bin == b, pct_lava := sum(sub[,prob_annot.lava])/nrow(sub) ]
# 	d[bin == b, pct_global := sum(sub[,prob_annot.global])/nrow(sub) ]
# 	print(d)
# }
# d

# x <- melt(d, id='bin')

# jpeg(paste0("/users/nicolerg/gtex-admix/plots/",tissue,"-pval-diff-log-v-annot-binary.jpg"),width=1600,height=1000,res=200)
# ggplot(x, aes(x=bin,y=value,fill=factor(variable),colour=factor(variable))) + 
# 	geom_bar(stat='identity',position='dodge') +
# 	labs(x='-log10(lava_pval) - -log10(global_pval)\nincreasing significance with LAVA ------->', y='Fraction of eQTLs with annotated top SNP',title=tissue) +
# 	scale_fill_manual(values=c(pct_global='blue',pct_lava='red'), 
# 		labels=c(pct_global='Global',pct_lava='LAVA'),name='Method') +
# 	scale_colour_manual(guide='none',
# 		values=c(pct_global='blue',pct_lava='red')) +
# 	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
# 	scale_x_discrete(breaks=bin_lab[,bins], labels=bin_lab[,labels]) +
# 	coord_cartesian(ylim=c(0,max(x[,value])+0.1)) +
# 	annotate('text',x=bin_lab[,bins],y=max(x[,value])+0.1,label=bin_lab[,counts])
# dev.off()

#########################################################################################

# # plot annotation of non-overlapping lead SNPs - only high significance 

# summ <- fread(infile, sep='\t', header=TRUE)
 
# # apply cuttoff
# summ <- summ[-log10(global_pval) > 6 & -log10(lava_pval) > 6]

# summ[,pval_diff := -log10(lava_pval) - -log10(global_pval)]

# # binarize 
# summ[,prob_annot.global := ifelse(global_annot > 0, 1, 0)]
# summ[,prob_annot.lava := ifelse(lava_annot > 0, 1, 0)]
# summ

# p_global <- sum(summ[,prob_annot.global]/nrow(summ))
# p_lava <- sum(summ[,prob_annot.lava]/nrow(summ))

# p_global
# p_lava 

# df <- data.frame(method=c('Global','LAVA'),p=c(p_global,p_lava))

# jpeg(paste0("/users/nicolerg/gtex-admix/plots/",tissue,"-annot-binary-highsig6.jpg"),width=400,height=600,res=150)
# ggplot(df, aes(x=method,y=p,fill=factor(method))) + 
# 	geom_bar(stat='identity',position='dodge') +
# 	labs(y='Fraction of eQTLs with annotated top SNP\n(p-val < 10e-5)',title=paste0(tissue,'\n(',nrow(summ),')')) +
# 	scale_fill_manual(values=c(Global='blue',LAVA='red'), guide=FALSE) +
# 	theme(axis.title.x=element_blank())
# dev.off()

#########################################################################################

# # make manhattan plots of the bins where gene counts are low using filtered allpairs file in conjunction with annotated top SNPs

# # get genes where counts are low

# get <- as.numeric(bin_lab[counts < 8, bins])
# select_genes <- summ[ bin %in% get, gene_id]

# # read in full annotations
# full <- fread(full, sep='\t', header=TRUE)
# full <- full[ gene_id %in% select_genes , ]

# ###############
# allpairs <- fread(cmd=sprintf("zcat %s",allp), sep='\t', header=TRUE)
# head(allpairs)
# allpairs <- allpairs[ gene_id %in% select_genes ]

# f <- function(x){
# 	return(as.numeric(unname(unlist(strsplit(x,'_')))[2]))
# }

# allpairs[,position := sapply(variant_id,f)]

# # save(allpairs,file="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/allp-filt.R")
# # ################
# # load(file="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/allp-filt.R")
# # merge full and allpairs 

# full <- full[,.(gene_id,variant_id,pval_nominal,method,annotation)]
# nrow(full)
# length(unique(full[,gene_id]))

# allp_global <- allpairs[,.(gene_id,variant_id,global_pval_nominal,position,tss_distance)]
# colnames(allp_global) <- c('gene_id','variant_id','pval_nominal','position','tss_distance')
# allp_global[,method:='global']

# allp_lava <- allpairs[,.(gene_id,variant_id,lava_pval_nominal,position,tss_distance)]
# colnames(allp_lava) <- c('gene_id','variant_id','pval_nominal','position','tss_distance')
# allp_lava[,method:='LAVA']

# allp <- data.table(rbind(allp_global,allp_lava))

# allp[,pval_nominal := signif(as.numeric(pval_nominal),8)]
# full[,pval_nominal := signif(as.numeric(pval_nominal),8)]

# merged <- merge(allp,full,by=c('gene_id','variant_id','method','pval_nominal'),all=TRUE)
# nrow(merged)

# merged[is.na(annotation), annotation := -1]

# ngenes <- length(select_genes)
# height <- 300*ngenes

# dtlist <- list()
# i <- 1
# # remove outliers
# for (g in unique(merged[,gene_id])){
# 	max <- max(merged[gene_id==g & annotation != -1, position])
# 	min <- min(merged[gene_id==g & annotation != -1, position])
# 	m <- mean(merged[gene_id==g,position])
# 	s <- sd(merged[gene_id==g,position])
# 	upper <- max(max, m+2*s)
# 	lower <- min(min, m-2*s)
# 	dtlist[[i]] <- merged[gene_id == g & position >= lower & position <= upper]
# 	i <- i + 1
# }

# m <- data.table(rbindlist(dtlist))
# paste(tissue,nrow(merged),nrow(m))

# jpeg(paste0('~/gtex-admix/plots/',tissue,'-top-snps-manhattan-nol-log.jpg'),width=1500,height=height,res=200)
# ggplot(m,aes(y=-log10(pval_nominal),x=position,colour=factor(method),shape=factor(annotation),alpha=factor(annotation))) +
# 	geom_point(fill='yellow') + 
# 	scale_colour_manual(values=c(LAVA='red',global='blue'),
# 		labels=c(global='Global',LAVA='LAVA'),
# 		name='Method') +
# 	scale_shape_manual(values=c('-1'=20,'1'=24,'0'=17),labels=c('1'='Top SNP,\nannotated','-1'='Other (p-val < 0.05)','0'='Top SNP,\nnot annotated'),name='Annotation') +
# 	scale_alpha_manual(values=c('0'=1,'1'=1,'-1'=0.3),guide=FALSE) +
# 	facet_wrap(~gene_id,scales='free',ncol=1) +
# 	labs(x='Position',y='-log10(nominal p-value)',title=tissue)
# dev.off()

#########################################################################################

# plot TSS distances for matching eGenes with non-overlapping lead SNPs

# full <- fread(full, sep='\t', header=TRUE)
# head(full)
# allpairs <- fread(cmd=sprintf("zcat %s",allp), sep='\t', header=TRUE)
# allpairs <- allpairs[,.(gene_id,variant_id,tss_distance)]

# tss <- merge(full, allpairs, by=c('gene_id','variant_id'))
# head(tss)

# tss_min <- data.table(gene_id=unique(tss[,gene_id]),lava=10e7,global=10e7)

# for (g in unique(tss[,gene_id])){
# 	sub <- tss[gene_id==g]
# 	lava <- sub[method=='LAVA']
# 	global <- sub[method=='global']
# 	lava_min <- lava[abs(tss_distance) == min(abs(tss_distance)), tss_distance] 
# 	global_min <- global[abs(tss_distance) == min(abs(tss_distance)), tss_distance]  

# 	tss_min[gene_id == g, lava := lava_min]
# 	tss_min[gene_id == g, global := global_min]
# }

# table(abs(tss_min[,global]) > abs(tss_min[,lava]))

# tss_min <- melt(tss_min, id='gene_id')

# jpeg(paste0('~/gtex-admix/plots/',tissue,'-best-tss.jpg'),width=1200,height=800,res=200)
# ggplot(tss_min, aes(x=value, fill=factor(variable))) +
# 	geom_density(alpha=0.2) +
# 	scale_fill_manual(values=c(global='blue',lava='red'),name='Method',
# 		labels=c(global='Global',lava='LAVA')) +
# 	labs(title=paste(tissue, '(same eGenes; non-overlapping top SNPs)'),x='Minimum TSS distance in top SNPs',y='Density') +
# 	geom_vline(xintercept=0,linetype="dashed")
# dev.off()

#########################################################################################

# # annotate SNPs with CADD scores (same eGenes, non-overlapping eSNPs)

# full <- fread(full, sep='\t', header=TRUE)
# cadd <- fread(caddmap, sep='\t', header=FALSE)
# colnames(cadd) <- c('variant_id','raw_score')
# head(cadd)

# # remove scores for variants that are not SNPs
# alt_f <- function(x){
# 	return(unname(unlist(strsplit(x,'_')))[4])
# }
# cadd[,alt:=sapply(variant_id,alt_f)]
# cadd <- cadd[nchar(alt)==1]

# # merge tables to assign CADD scores to variants
# cadd[,alt:=NULL]
# m <- merge(full, cadd, by='variant_id',all.x=TRUE)

# summ <- fread(infile, sep='\t', header=TRUE)
# summ <- summ[,.(gene_id,global_pval,lava_pval)]
# summ[,pval_diff := global_pval - lava_pval]
# summ[,global_pval := NULL]
# summ[,lava_pval := NULL]

# master <- merge(m,summ,by='gene_id',all.x=TRUE)
# head(master)

# # divide data into bins

# bins <- bin(master[,pval_diff],nbins=10,method='length',labels=1:10)
# master[,bin := bins]
# head(master)

# bin_labels <- bin(master[,pval_diff],nbins=10,method='length')
# bin_lab <- data.table(bins=bins, labels=bin_labels)
# bin_lab <- bin_lab[!duplicated(bins)]
# bin_lab <- bin_lab[order(bins)]
# x <- data.table(table(master[,bin]))
# x <- x[N != 0]
# bin_lab[,counts := x[,N]] # V1, N
# print(bin_lab)

# y <- max(master[,raw_score],na.rm=TRUE)+0.1

# jpeg(paste0("/users/nicolerg/gtex-admix/plots/",tissue,"-pval-diff-v-cadd.jpg"),width=1600,height=1000,res=200)
# ggplot(master, aes(x=bin,y=raw_score,colour=factor(method))) + 
# 	geom_boxplot() +
# 	labs(x='Difference in nominal p-value of top SNPs (global - LAVA)', y='Raw CADD scores of tied top SNPs',title=tissue) +
# 	scale_colour_manual(values=c(global='blue',LAVA='red'),name='Method') +
# 	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
# 	scale_x_discrete(breaks=bin_lab[,bins], labels=bin_lab[,labels]) +
# 	annotate('text',x=bin_lab[,bins],y=y,label=bin_lab[,counts])
# dev.off()


#########################################################################################

# histogram of p-values from filtered allpairs

allpairs <- fread(cmd=sprintf("zcat %s",allp), sep='\t', header=TRUE)
allp <- allpairs[,.(global_pval_nominal,lava_pval_nominal)]
allp[,global_pval_nominal := -log10(global_pval_nominal)]
allp[,lava_pval_nominal := -log10(lava_pval_nominal)]

melt_allp <- melt(allp,measure.vars=c('global_pval_nominal','lava_pval_nominal'))

jpeg(paste0("/users/nicolerg/gtex-admix/plots/",tissue,"-pval-hist.jpg"),width=1600,height=1000,res=200)
ggplot(melt_allp, aes(x=value,fill=factor(variable))) + 
	geom_histogram(bins=100,alpha=0.5,position='identity') +
	labs(x='-log10(p-value)', y='Count',title=tissue) +
	scale_fill_manual(values=c(global_pval_nominal='blue',lava_pval_nominal='red'),
		labels=c(global_pval_nominal='Global',lava_pval_nominal='LAVA'),name='Method') 
dev.off()


#########################################################################################

# q-q plot of LAVA vs global p-values (p < 0.05)

jpeg(paste0("/users/nicolerg/gtex-admix/plots/",tissue,"-qq.jpg"),width=1000,height=1000,res=200)
qqplot(x=allp[,global_pval_nominal], y=allp[,lava_pval_nominal], 
	main=tissue, xlab='-log10(global p-value)', ylab='-log10(LAVA p-value)', pch=19, cex=0.5)
abline(a=0,b=1,lty='dotted')
dev.off()

#########################################################################################



# # annotate SNPs with allele-specific binding (same eGenes, non-overlapping eSNPs)

# hg19 <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/ng.3432-S7'

