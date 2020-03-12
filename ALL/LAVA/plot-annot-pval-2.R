library(ggplot2)
library(data.table)
require('OneR')
library(reshape)
library(stats)

args <- commandArgs(trailingOnly=TRUE)
tissue <- args[1]
infile <- args[2]
allp <- args[3]

# infile <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/Lung.diff.SNP.annot.tsv'
# tissue <- 'Lung'
# full <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/Lung.diff.eSNP.annotated.tsv'
# allp <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/Lung-LAVA-global-merged-allpairs-filt.tsv.gz'

#########################################################################################

# plot annotation of different eGenes 

annot <- fread(infile,sep='\t',header=TRUE)

dt <- data.table(gene_id=unique(annot[,gene_id]),method='None',count=0,annotation=0)
i <- 1
# collapse on gene
for (g in unique(annot[,gene_id])){
	sub <- annot[gene_id == g]
	m <- unique(sub[,method])
	total <- nrow(sub)
	sum <- sum(sub[,annotation])

	dt[gene_id == g, method := m]
	dt[gene_id == g, count := total]
	dt[gene_id == g, annotation := ifelse(sum > 0,1,0)]
}

head(dt)
head(dt)


g_total <- nrow(dt[method=='global',])
g_annot <- sum(dt[method=='global',annotation])

l_total <- nrow(dt[method=='LAVA'])
l_annot <- sum(dt[method=='LAVA',annotation])

df <- data.table(method=c('global','LAVA'),
	total_egenes=c(g_total,l_total),
	total_annot_egenes=c(g_annot,l_annot))

df

df[,total_not_annot := total_egenes - total_annot_egenes]

df <- melt(df, id=c('method','total_egenes'))

jpeg(paste0("/users/nicolerg/gtex-admix/plots/",tissue,"-diff-egene-annot.jpg"),width=1600,height=1000,res=200)
ggplot(df, aes(x=factor(method),y=value,fill=factor(variable))) + 
	geom_bar(stat='identity') +
	scale_fill_manual(values=c(total_annot_egenes='green',total_not_annot='gray'),
		labels=c(total_annot_egenes='eGenes with annotated top SNP',
			total_not_annot='eGenes with no annotated top SNP'),name='') +
	labs(y='N unique eGenes',x='Method',title=tissue)
dev.off()









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
