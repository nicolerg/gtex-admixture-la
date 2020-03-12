library(data.table)
library(ggplot2)

indir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/"

# plot effect size in GEUVADIS eQTLs ========================================================================================================================

eurgeu <- "EUR373.gene.cis.FDR5.best.rs137.hg38.txt.gz"
afrgeu <- "YRI89.gene.cis.FDR5.best.rs137.hg38.txt.gz"

eur <- fread(input=sprintf("zcat %s",eurgeu), sep='\t', header=TRUE)
afr <- fread(input=sprintf("zcat %s",afrgeu), sep='\t', header=TRUE)

merged <- merge(eur, afr, by="GENE_ID")
print(nrow(merged))

# color eQTL genes with significantly different effect sizes 
mean <- mean(abs(c(merged[,rvalue.x], merged[,rvalue.y])))
merged[, label := as.integer(abs(rvalue.x - rvalue.y) >= mean) ]

pdf("~/gtex-admix/plots/geu-eqtl-effect-sizes.pdf", width=6, height=6)
ggplot(merged, aes(x=rvalue.x,y=rvalue.y,colour=factor(label))) +
	geom_point() +
	labs(x="EUR eQTL r-value",y="AFR eQTL r-value") +
	theme_bw() +
	scale_colour_manual(values=c("1"="red","0"="black")) +
	theme(legend.position="none")
dev.off()

# write out list of eQTLs with higher effect size in AFR or EUR

merged <- merged[label==1]

colnames(merged) <- gsub(".x",".EUR",colnames(merged),fixed=TRUE)
colnames(merged) <- gsub(".y",".AFR",colnames(merged),fixed=TRUE)

highafr <- merged[rvalue.AFR > rvalue.EUR]
higheur <- merged[rvalue.EUR > rvalue.AFR]

write.table(highafr,paste0(indir,'GEU-eQTLS-highAFR-effect.tsv'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(higheur,paste0(indir,'GEU-eQTLS-highEUR-effect.tsv'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
