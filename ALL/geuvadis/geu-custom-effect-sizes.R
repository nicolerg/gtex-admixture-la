library(data.table)
library(ggplot2)
library(qvalue)

# indir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/"

# # pval correction
# # merge AFR allpairs with EUR eQTLs, AFR eQTLs with EUR allpairs

# afrall <- paste0(indir, 'AFR.geu.allpairs.txt.gz')
# eurall <- paste0(indir, 'EUR.geu.allpairs.txt.gz')
# afrg <- paste0(indir, 'AFR.geu.permuted.txt.gz')
# eurg <- paste0(indir, 'EUR.geu.permuted.txt.gz')

# # ID egenes ===================================================================================================

# eurg <- fread(input=sprintf("zcat %s", eurg), sep=' ', header=TRUE)
# afrg <- fread(input=sprintf("zcat %s", afrg), sep=' ', header=TRUE)

# eurg[,qval:=qvalue(pval_beta)$qvalues]
# afrg[,qval:=qvalue(pval_beta)$qvalues]

# eurg <- eurg[qval < 0.05]
# print(paste("N EUR eQTLs:",nrow(eurg)))
# afrg <- afrg[qval < 0.05]
# print(paste("N AFR eQTLs:",nrow(afrg)))

# # get matching pairs ===========================================================================================

# eurall <- fread(input=sprintf("zcat %s",eurall), sep=' ', header=TRUE)
# afrall <- fread(input=sprintf("zcat %s",afrall), sep=' ', header=TRUE)

# eurbest <- merge(eurg, afrall, by=c("gene_id","variant"), all.x=TRUE, suffixes=c(".EUR",".AFR"))
# afrbest <- merge(afrg, eurall, by=c("gene_id","variant"), all.x=TRUE, suffixes=c(".AFR",".EUR"))

# # write out different gene sets

# # # set 1: AFR eQTLs
# # afreqtl <- afrg[,gene_id]
# # afreqtl <- gsub("\\..*","",afreqtl)
# # write.table(afreqtl,paste0(indir,"afr.eqtl.genes.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # # set 2: EUR eQTLs
# # eureqtl <- eurg[,gene_id]
# # eureqtl <- gsub("\\..*","",eureqtl)
# # write.table(eureqtl,paste0(indir,"eur.eqtl.genes.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # # set 3: AFR eQTLs with non-zero EUR effect size
# # afrfilt <- afrbest[slope.EUR != 0, gene_id]
# # afrfilt <- gsub("\\..*","",afrfilt)
# # write.table(afrfilt,paste0(indir,"afr.eqtl.genes.filt.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # # set 4: EUR eQTLs with non-zero AFR effect size
# # eurfilt <- eurbest[slope.AFR != 0, gene_id]
# # eurfilt <- gsub("\\..*","",eurfilt)
# # write.table(eurfilt,paste0(indir,"eur.eqtl.genes.filt.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# afrbest[,colour:="AFR"]
# eurbest[,colour:="EUR"]
# allegene <- data.table(rbind(eurbest, afrbest))

# head(allegene)
# save(allegene, file=paste0(indir, "AFR-EUR-eqtls-merged.RData"))

# plot effect size in GEUVADIS eQTLs ===========================================================================

# load(file=paste0(indir, "AFR-EUR-eqtls-merged.RData"))

# allegene[, sigpoint := 0 ]
# allegene[ abs(slope.AFR) > (abs(slope.EUR) + sqrt(abs(slope.EUR))), sigpoint := 1 ]
# allegene[ abs(slope.AFR) < (abs(slope.EUR) - sqrt(abs(slope.EUR))), sigpoint := 1 ]

# allegene[ slope.AFR > exp(slope.EUR) - 0.5 & slope.AFR > exp(-slope.EUR) - 0.5 , sigpoint := 1]
# allegene[ slope.AFR < -exp(slope.EUR) + 0.5 & slope.AFR < -exp(-slope.EUR) + 0.5 , sigpoint := 1]
# allegene[ slope.AFR < log(slope.EUR+0.5) & slope.AFR > -log(slope.EUR+0.5) , sigpoint := 1]
# allegene[ slope.AFR < log(-slope.EUR+0.5) & slope.AFR > -log(-slope.EUR+0.5) , sigpoint := 1]

# sig <- allegene[ sigpoint == 1 ]
# sig[,sigpoint:=NULL]

# write out different gene sets

# # set 1: AFR eQTLs
# afreqtl <- sig[colour=="AFR",gene_id]
# afreqtl <- gsub("\\..*","",afreqtl)
# write.table(afreqtl,paste0(indir,"afr.eqtl.genes.diff.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # set 2: EUR eQTLs
# eureqtl <- sig[colour=="EUR",gene_id]
# eureqtl <- gsub("\\..*","",eureqtl)
# write.table(eureqtl,paste0(indir,"eur.eqtl.genes.diff.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # set 3: AFR eQTLs with non-zero EUR effect size
# afrfilt <- sig[colour=="AFR" & slope.EUR != 0,gene_id]
# afrfilt <- gsub("\\..*","",afrfilt)
# write.table(afrfilt,paste0(indir,"afr.eqtl.genes.filt.diff.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # set 4: EUR eQTLs with non-zero AFR effect size
# eurfilt <- sig[colour=="EUR" & slope.AFR != 0, gene_id]
# eurfilt <- gsub("\\..*","",eurfilt)
# write.table(eurfilt,paste0(indir,"eur.eqtl.genes.filt.diff.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # EUR eQTLs with 0 slope in AFR
# afrmono <- allegene[colour=="EUR" & slope.AFR == 0]
# afrmono[,sigpoint := NULL]
# write.table(afrmono,paste0(indir,"eur.eqtl.afr.0slope.txt"),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)

# regions <- data.frame(x1 = seq(0,max(allegene[,slope.EUR]),length.out=50),
# 	x2 = seq(0.5,max(allegene[,slope.EUR]),length.out=50))
# regions$y11 <- exp(regions$x1) - 0.5
# regions$y12 <- -exp(regions$x1) + 0.5
# regions$y31 <- log(regions$x2+0.5)

# jpeg("~/gtex-admix/plots/geu-custom-eqtl-effect-sizes.jpg", width=1200, height=1000, res=200)
# ggplot() +
# 	geom_point(data=allegene,aes(x=slope.EUR,y=slope.AFR,colour=factor(colour)),alpha = 0.7, size = 0.5) +
# 	labs(x="EUR effect size",y="AFR effect size") +
# 	scale_colour_manual(values=c("EUR"="orange","AFR"="blue"),name="eQTL population") +
# 	guides(colour = guide_legend(override.aes = list(size=2))) +
# 	geom_point(data=sig, aes(x=slope.EUR,y=slope.AFR,fill=factor(colour)), size=1, alpha=1, colour="black", shape=21) +
# 	scale_fill_manual(name="Different effect sizes",values=c("EUR"="orange","AFR"="blue"),guide='none') +
# 	geom_line(data=regions,aes(x=x1,y=y11),linetype="dashed") +
# 	geom_line(data=regions,aes(x=-x1,y=y11),linetype="dashed") +
# 	geom_line(data=regions,aes(x=x1,y=y12),linetype="dashed") +
# 	geom_line(data=regions,aes(x=-x1,y=y12),linetype="dashed") +
# 	geom_line(data=regions,aes(x=x2,y=y31),linetype="dashed") +
# 	geom_line(data=regions,aes(x=x2,y=-y31),linetype="dashed") +
# 	geom_line(data=regions,aes(x=-x2,y=y31),linetype="dashed") +
# 	geom_line(data=regions,aes(x=-x2,y=-y31),linetype="dashed") 
# dev.off()

# for EUR into AFR 

# args <- commandArgs(trailingOnly=TRUE)
# indir <- paste0(args[1],'/')
# # indir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/"

# afrall <- paste0(indir, 'AFR.geu.allpairs.txt.gz')
# eurall <- paste0(indir, 'EUR.geu.allpairs.txt.gz')
# afrg <- paste0(indir, 'AFR.geu.permuted.txt.gz')
# eurg <- paste0(indir, 'EUR.geu.permuted.txt.gz')

# # ID egenes ===================================================================================================

# eurg <- fread(input=sprintf("zcat %s", eurg), sep=' ', header=TRUE)
# afrg <- fread(input=sprintf("zcat %s", afrg), sep=' ', header=TRUE)

# eurg[,qval:=qvalue(pval_beta)$qvalues]
# afrg[,qval:=qvalue(pval_beta)$qvalues]

# eurg <- eurg[qval < 0.05]
# print(paste("N EUR eQTLs:",nrow(eurg)))
# afrg <- afrg[qval < 0.05]
# print(paste("N AFR eQTLs:",nrow(afrg)))

# # get matching pairs ===========================================================================================

# eurall <- fread(input=sprintf("zcat %s",eurall), sep=' ', header=TRUE)
# afrall <- fread(input=sprintf("zcat %s",afrall), sep=' ', header=TRUE)

# eurbest <- merge(eurg, afrall, by=c("gene_id","variant"), all.x=TRUE, suffixes=c(".EUR",".AFR"))
# afrbest <- merge(afrg, eurall, by=c("gene_id","variant"), all.x=TRUE, suffixes=c(".AFR",".EUR"))

# # write out different gene sets

# # # set 1: AFR eQTLs
# # afreqtl <- afrg[,gene_id]
# # afreqtl <- gsub("\\..*","",afreqtl)
# # write.table(afreqtl,paste0(indir,"afr.eqtl.genes.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # # set 2: EUR eQTLs
# # eureqtl <- eurg[,gene_id]
# # eureqtl <- gsub("\\..*","",eureqtl)
# # write.table(eureqtl,paste0(indir,"eur.eqtl.genes.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # # set 3: AFR eQTLs with non-zero EUR effect size
# # afrfilt <- afrbest[slope.EUR != 0, gene_id]
# # afrfilt <- gsub("\\..*","",afrfilt)
# # write.table(afrfilt,paste0(indir,"afr.eqtl.genes.filt.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # # set 4: EUR eQTLs with non-zero AFR effect size
# # eurfilt <- eurbest[slope.AFR != 0, gene_id]
# # eurfilt <- gsub("\\..*","",eurfilt)
# # write.table(eurfilt,paste0(indir,"eur.eqtl.genes.filt.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# # get genes with 0 AFR slope
# afr0 <- eurbest[slope.AFR == 0]
# if(nrow(afr0) > 0){
# 	write.table(afr0,paste0(indir,"EUR.eqtl.AFR.0slope.txt"),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
# }

# # get genes with 0 EUR slope
# eur0 <- afrbest[slope.EUR == 0]
# if(nrow(eur0) > 0){
# 	write.table(eur0,paste0(indir,"AFR.eqtl.EUR.0slope.txt"),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
# }

# afrbest[,colour:="AFR"]
# eurbest[,colour:="EUR"]
# allegene <- data.table(rbind(eurbest, afrbest))

# head(allegene)
# save(allegene, file=paste0(indir, "AFR-EUR-eqtls-merged.RData"))

# load(file=paste0(indir, "AFR-EUR-eqtls-merged.RData"))
# allegene[, sigpoint := 0 ]
# allegene[,diff:=abs(slope.EUR - slope.AFR)]
# allegene <- allegene[order(diff,decreasing=TRUE)]
# #n <- round(0.05*nrow(allegene),0)
# allegene[1:100,sigpoint := 1]

# if(grepl("EUR-into-AFR",indir)){
# 	outfile <- '~/gtex-admix/plots/figures/geu-eur-into-afr.pdf'
# 	lims <- c(min(allegene[,slope.AFR]),max(allegene[,slope.AFR]))
# } else {
# 	outfile <- '~/gtex-admix/plots/figures/geu-afr-into-eur.pdf'
# 	lims <- c(min(allegene[,slope.EUR]),max(allegene[,slope.EUR]))
# }

# pdf(outfile,width=4,height=4)
# ggplot(data=NULL) +
# 	geom_point(data=allegene,aes(x=slope.EUR,y=slope.AFR,colour=factor(colour)),alpha=0.7,size=0.5) +
# 	labs(x="EUR effect size",y="AFR effect size") +
# 	scale_colour_manual(values=c("EUR"="#0000FF","AFR"="#FF9900"),name="eQTL population") +
# 	guides(colour = guide_legend(override.aes = list(size=2))) +
# 	geom_point(data=allegene[sigpoint==1], aes(x=slope.EUR,y=slope.AFR,fill=factor(colour)), size=1, alpha=1, colour="black", shape=21) +
# 	scale_fill_manual(name="Different effect sizes",values=c("EUR"="#0000FF","AFR"="#FF9900"),guide='none') +
# 	ylim(lims) +
# 	xlim(lims) +
# 	theme_bw() +
# 	theme(legend.position=c(0.78,0.15),legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"))
# dev.off()

# now do it for AFR into EUR 

# let's simplify this ====================================================================================================

mergeres <- function(qtl,all,suff){
	qtl <- fread(qtl, sep='\t',header=TRUE)
	all <- fread(input=sprintf("zcat %s",all),sep=' ',header=TRUE)
	merged <- merge(qtl,all,by=c('gene_id','variant'),suffixes=suff)
	return(merged)
}

# EUR into AFR
indir <- '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR-into-AFR/'
qtl <- paste0(indir,'EUR.geu.egenes.5fdr.txt')
all <- paste0(indir,'AFR.geu.allpairs.txt.gz')
eur2afr <- mergeres(qtl,all,c(".EUR",".AFR"))

# AFR into EUR
indir <- '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR/'
qtl <- paste0(indir,'AFR.geu.egenes.5fdr.txt')
all <- paste0(indir,'EUR.geu.allpairs.txt.gz')
afr2eur <- mergeres(qtl,all,c(".AFR",".EUR"))

eur2afr[,pop:='EUR2AFR']
afr2eur[,pop:='AFR2EUR']

m <- rbind(eur2afr,afr2eur)

m[, sigpoint := 0 ]
m[,diff:=abs(slope.EUR - slope.AFR)]
m <- m[order(diff,decreasing=TRUE)]
m[1:100,sigpoint := 1]

write.table(m,'/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/merged-eqtls.tsv',sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
# m <- fread("/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/merged-eqtls.tsv", sep='\t', header=TRUE)

# pdf("~/gtex-admix/plots/figures/geu-eqtl-by-pop.pdf",width=4,height=4)
# ggplot(data=NULL) +
# 	geom_point(data=m,aes(x=slope.EUR,y=slope.AFR,colour=factor(pop)),alpha=0.7,size=0.5) +
# 	labs(x="EUR effect size",y="AFR effect size") +
# 	scale_colour_manual(values=c(EUR2AFR="#0000FF",AFR2EUR="#FF9900"),name="eQTL population",
# 		labels=c(EUR2AFR='EUR into AFR', AFR2EUR='AFR into EUR'),guide='legend') +
# 	geom_point(data=m[sigpoint==1], aes(x=slope.EUR,y=slope.AFR,fill=factor(pop),shape=factor(sigpoint)), size=1, alpha=1, colour="black") +
# 	scale_shape_manual(values=21,labels='Top 100 most\ndifferent') +
# 	scale_fill_manual(name="Different effect sizes",values=c(EUR2AFR="#0000FF",AFR2EUR="#FF9900"),guide='none') +
# 	guides(colour = guide_legend(override.aes = list(size=2),order=1,keywidth=unit(0.9,"cm")),
# 		shape = guide_legend(override.aes = list(size=2),order=0,title=NULL,label.vjust=1.5,keywidth=unit(1,"cm"))) +
# 	# ylim(lims) +
# 	# xlim(lims) +
# 	theme_bw() +
# 	theme(legend.position=c(0.79,0.23),legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"))
# dev.off()

