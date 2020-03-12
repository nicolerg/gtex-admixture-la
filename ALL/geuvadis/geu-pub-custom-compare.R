# compare GEUVADIS r-values to custom eQTL calling betas 

library(data.table)
library(ggplot2)

# first with AFR
chr1allp <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR.chr1.geu.allpairs.txt.gz"
geupub <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/YRI89.gene.cis.FDR5.all.rs137.txt.gz"

allp <- fread(input=sprintf("zcat %s",chr1allp), sep=' ', header=TRUE)
allp <- allp[pval_nominal < 0.05]

# gene_id variant tss_distance pval_nominal slope
# ENSG00000223972.4 rs367896724 -1692 0.53841 0.154343
geup <- fread(input=sprintf("zcat %s",geupub), sep='\t', header=TRUE)
# SNP_ID ID GENE_ID PROBE_ID CHR_SNP CHR_GENE SNPpos TSSpos distance rvalue pvalue log10pvalue

geu <- geup[,list(GENE_ID,SNP_ID,CHR_GENE,rvalue)]
geu <- geu[CHR_GENE == 1]

geugenes <- unique(geu[,GENE_ID])
mygenes <- unique(allp[,gene_id])

print(length(geugenes))
print(length(mygenes))

# merge results
interx <- merge(allp, geu, by.x=c("gene_id","variant"), by.y=c("GENE_ID","SNP_ID"))
head(interx)

print('N my pairs:',quote=FALSE)
nrow(allp)
print('N GEU pairs:',quote=FALSE)
nrow(geu)
print('N pairs in intersection:', quote=FALSE)
nrow(interx)

# jpeg("~/gtex-admix/plots/geuvadis-effect-replication-afr.jpg",width=1200,height=1200,res=200)
# ggplot(interx, aes(x=slope, y=rvalue)) +
# 	geom_point(size=0.5) +
# 	theme_bw() +
# 	labs(x="beta (custom eQTL calling pipeline)",y="rvalue (published data)")
# dev.off()

# now with EUR

chr1allp <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR.chr1.geu.allpairs.txt.gz"
geupub <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/EUR373.gene.cis.FDR5.all.rs137.txt.gz"

allpeur <- fread(input=sprintf("zcat %s",chr1allp), sep=' ', header=TRUE)
allpeur <- allpeur[pval_nominal < 0.05]

# gene_id variant tss_distance pval_nominal slope
# ENSG00000223972.4 rs367896724 -1692 0.53841 0.154343
geup <- fread(input=sprintf("zcat %s",geupub), sep='\t', header=TRUE)
# SNP_ID ID GENE_ID PROBE_ID CHR_SNP CHR_GENE SNPpos TSSpos distance rvalue pvalue log10pvalue

geueur <- geup[,list(GENE_ID,SNP_ID,CHR_GENE,rvalue)]
geueur <- geueur[CHR_GENE == 1]

geugenes <- unique(geueur[,GENE_ID])
mygenes <- unique(allpeur[,gene_id])

print(length(geugenes))
print(length(mygenes))

# merge results
xeur <- merge(allpeur, geueur, by.x=c("gene_id","variant"), by.y=c("GENE_ID","SNP_ID"))
head(xeur)

print('N my pairs:',quote=FALSE)
nrow(allpeur)
print('N GEU pairs:',quote=FALSE)
nrow(geueur)
print('N pairs in intersection:', quote=FALSE)
nrow(xeur)

# jpeg("~/gtex-admix/plots/geuvadis-effect-replication-eur.jpg",width=1200,height=1200,res=200)
# ggplot(xeur, aes(x=slope, y=rvalue)) +
# 	geom_point(size=0.5) +
# 	theme_bw() +
# 	labs(x="beta (custom eQTL calling pipeline)",y="rvalue (published data)")
# dev.off()

xeur[,pop:="EUR"]
interx[,pop:="AFR"]

master <- rbind(xeur, interx)

jpeg("~/gtex-admix/plots/geuvadis-effect-replication-all.jpg",width=1000,height=800,res=200)
ggplot(master, aes(x=slope, y=rvalue, colour=pop)) +
	geom_point(size=0.6,alpha=0.5) +
	theme_bw() +
	labs(x="beta (custom eQTL calling pipeline)",y="rvalue (published data)") + 
	scale_colour_manual(name="Population",values=c(EUR="blue",AFR="orange"))
dev.off()

# # compare effect sizes from different populations ==============================================================

chr1eur <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR.chr1.geu.allpairs.txt.gz"
allpeur <- fread(input=sprintf("zcat %s",chr1eur), sep=' ', header=TRUE)

cutoff <- 0.05 / length(unique(allpeur[,gene_id]))
eur <- allpeur[pval_nominal < 0.95]

eur[,sig_eur:=ifelse(pval_nominal < cutoff,1,0)]

notsig <- eur[sig_eur==0]
sig <- eur[sig_eur==1]

# pval vs effect size in EUR only
jpeg("~/gtex-admix/plots/geuvadis-eur-pval-effect.jpg",width=800,height=800,res=200)
ggplot(data=NULL, aes(y=-log(pval_nominal), x=slope)) +
	geom_point(data=notsig,size=0.6,alpha=0.1,colour="black") +
	geom_point(data=sig,size=0.6,alpha=0.1,colour="red") +
	theme_bw() +
	labs(y="-log(EUR nominal pval)",x="EUR effect size")
dev.off()

chr1afr <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR.chr1.geu.allpairs.txt.gz"
allpafr <- fread(input=sprintf("zcat %s",chr1afr), sep=' ', header=TRUE)

mergeallp <- merge(allpafr, allpeur, by=c("gene_id","variant","tss_distance"), suffixes = c(".AFR",".EUR"))

head(mergeallp)
write.table(mergeallp,file="/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/merged-all-pairs.tsv",sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)

mergeallp <- fread("/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/merged-all-pairs.tsv",sep='\t',header=TRUE)

nrow(mergeallp)
cutoff <- 0.05 / length(unique(mergeallp[,gene_id]))
print(cutoff)
mergeallp <- mergeallp[pval_nominal.EUR < 1.0 & pval_nominal.AFR < 1.0]
nrow(mergeallp)

mergeallp[,sig_eur:=ifelse(pval_nominal.EUR < cutoff,1,0)]
head(mergeallp)

notsig <- mergeallp[sig_eur==0]
sig <- mergeallp[sig_eur==1]

jpeg("~/gtex-admix/plots/geuvadis-custom-pop-effect-sig.jpg",width=800,height=800,res=200)
ggplot(data=NULL, aes(x=slope.EUR, y=slope.AFR)) +
	geom_point(data=notsig,size=0.6,alpha=0.1,colour="black") +
	geom_point(data=sig,size=0.6,alpha=0.1,colour="red") +
	theme_bw() +
	labs(x="EUR gene-SNP pair effect size",y="AFR gene-SNP pair effect size") +
	theme(legend.position="none")
dev.off()

q()

# see what it looks like in GTEx =================================================================================

inf <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/Cells_Cultured_fibroblasts.allpairs.txt.gz"
lcl <- fread(input=sprintf("zcat %s",inf),sep='\t',header=TRUE)

lcl <- lcl[grepl("chr1_",variant_id)]
cutoff <- 0.05 / (length(unique(lcl[,gene_id])))
print(cutoff)
print(nrow(lcl))
lcl <- lcl[pval_nominal < 0.95]
lcl[,sig:=ifelse(pval_nominal < cutoff,1,0)]
print(nrow(lcl))

notsig <- lcl[sig==0]
sig <- lcl[sig==1]

jpeg("~/gtex-admix/plots/gtex-lcl-pval-effect.jpg",width=800,height=800,res=200)
ggplot(data=NULL, aes(y=-log(pval_nominal), x=slope)) +
	geom_point(data=notsig,size=0.6,alpha=0.1,colour="black") +
	geom_point(data=sig,size=0.6,alpha=0.1,colour="red") +
	theme_bw() +
	labs(y="-log(nominal pval)",x="effect size")
dev.off()