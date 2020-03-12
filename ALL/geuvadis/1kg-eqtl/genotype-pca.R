library(SNPRelate)
library(data.table)
library(ggplot2)
library(plyr)

pcadir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/pca/"
vcfdir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/genotype/"
peerdir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/peer/"
covdir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/cov/"
tg.pop.codes <- "~/gtex-admix/metadata/geuvadis-cov.tsv"
meta <- read.delim(tg.pop.codes, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# # identify outlier
# load(file=paste0(pcadir,"tg.AFR.pca.RData"))

# df <- data.frame(sample.id = intersect_pca$sample.id,
# 	 Eig1 = intersect_pca$eigenvect[,1])
# df$sample.id[df$Eig1 > 0.6]
# # outlier: NA19118
# q()

for(pop in c("EUR","AFR")){
	# vcfin <- paste0(pcadir,"tg.",pop,".subset.recode.vcf.gz")
	# gds <- paste0(pcadir,"tg.",pop,"subset.gds")
	# snpgdsVCF2GDS(vcfin, gds,  method="biallelic.only")

	# # perform PCA
	# genofile <- snpgdsOpen(gds)
	# system("echo 'Performing PCA...'")
	# intersect_pca <- snpgdsPCA(genofile, num.thread=12)
	# save(intersect_pca, file=paste0(pcadir,"tg.",pop,".pca.RData"))
	load(file=paste0(pcadir,"tg.",pop,".pca.RData"))

	vcfdf <- data.frame(sample.id = intersect_pca$sample.id,
    Eig1 = intersect_pca$eigenvect[,1],
    Eig2 = intersect_pca$eigenvect[,2],
    Eig3 = intersect_pca$eigenvect[,3],
    stringsAsFactors = FALSE)

    vcfdf$pop <- meta$Population[match(vcfdf$sample.id, meta$Sample)]
    obj <- paste0(pop,"_pcdf")
    assign(obj, vcfdf)
}

# group.colors <- c(CEU="#0000FF",
# 	FIN="#00C5CD",
# 	GBR="#00EBFF",
# 	TSI="#00008B")

# pdf("~/gtex-admix/plots/geuvadis/EUR-eig12.pdf", height=6, width=8)
# ggplot(EUR_pcdf, aes(x=Eig1, y=Eig2, colour=pop)) +
# 	geom_point(alpha=1, size=1) + 
# 	theme(legend.position="right") +
# 	scale_colour_manual(values=group.colors,name="Euro Sub-Pop")
# dev.off()

# pdf("~/gtex-admix/plots/geuvadis/AFR-eig12.pdf", height=6, width=8)
# ggplot(AFR_pcdf, aes(x=Eig1, y=Eig2)) +
# 	geom_point(alpha=1, size=1, colour="#FFB933")
# dev.off()

for (pop in c("EUR", "AFR")){
	cov <- meta[meta$Super == pop,c("Sample","Gender")]
	if (pop == "EUR"){
		cov <- merge(cov, EUR_pcdf, by.x="Sample", by.y="sample.id")
		cov$Eig3 <- NULL
	} else if (pop == "AFR"){
		cov <- merge(cov, AFR_pcdf, by.x="Sample", by.y="sample.id")
	}
	cov$pop <- NULL
	# PEER factors
	peer <- read.delim(paste0(peerdir,"GV-",pop,".PEER_covariates.txt"), sep='\t', header=TRUE, stringsAsFactors=FALSE)
	ids <- peer$ID 
	peer$ID <- NULL
	subjs <- colnames(peer)
	peer <- data.frame(t(peer),stringsAsFactors=FALSE)
	colnames(peer) <- ids
	peer$SUBJECT <- subjs
	merged <- merge(cov, peer, by.x="Sample", by.y="SUBJECT")
	subjs <- merged$Sample
	merged$Sample <- NULL
	cols <- colnames(merged)
	merged <- data.frame(t(merged),stringsAsFactors=FALSE)
	colnames(merged) <- subjs
	merged$ID <- cols
	merged <- merged[,c("ID",subjs)]
	rownames(merged) <- NULL
	write.table(merged, paste0(covdir,pop,".covariates.tsv"), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
}
