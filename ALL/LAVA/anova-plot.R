library(data.table)
library(ggplot2)


outdir <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/geuvadis-curated"
anovapairs <- paste0(outdir,"/Artery_Tibial-curated-lava-anova-allpairs-chr22.merged.tsv.gz")
# anova_pval      anova_global_rss        anova_local_rss
allpairs <- fread(input=sprintf("zcat %s", anovapairs), sep='\t', header=TRUE)

allpairs[,pos := as.numeric(unlist(lapply(strsplit(variant_id,split='_'), '[', 2))) ]

# allpairs$pos <- as.numeric(unlist(lapply(strsplit(allpairs$variant_id,split='_'), '[', 2)))
# allpairs$anova_pval_adj <- -log10(p.adjust(allpairs$anova_pval, method="BH"))

# allpairs <- data.frame(allpairs, stringsAsFactors=FALSE)
allpairs <- na.omit(allpairs)
head(allpairs)

bonthresh <- -log10(0.05/nrow(allpairs))
print(bonthresh)

allpairs[,global_pval := -log10(as.numeric(global_pval))]
allpairs[,local_pval := -log10(as.numeric(local_pval))]
allpairs[,anova_pval := -log10(as.numeric(anova_pval))]
allpairs[,afr_pval := -log10(as.numeric(afr_pval))]

# allpairs <- allpairs[order(allpairs$anova_pval_adj, decreasing=TRUE),]
# sig <- unique(allpairs[1:10000,1])[1:10]
# sig

# sig <- c("ENSG00000223745.7",
# 	"ENSG00000281741.2",
# 	"ENSG00000134202.10",
# 	"ENSG00000272030.1")

# allpairs$g1 <- NA
# allpairs$g2 <- NA
# allpairs$g3 <- NA
# allpairs$g4 <- NA

# allpairs$g1[grep(sig[1], allpairs$gene_id)] <- allpairs$anova_pval[grep(sig[1], allpairs$gene_id)]
# allpairs$g2[grep(sig[2], allpairs$gene_id)] <- allpairs$anova_pval[grep(sig[2], allpairs$gene_id)]
# allpairs$g3[grep(sig[3], allpairs$gene_id)] <- allpairs$anova_pval[grep(sig[3], allpairs$gene_id)]
# allpairs$g4[grep(sig[4], allpairs$gene_id)] <- allpairs$anova_pval[grep(sig[4], allpairs$gene_id)]

head(allpairs)

pdf('~/gtex-admix/plots/Artery_Tibial-geuvadis-curated-chr22.merged.pdf', width=10, height=6)
g <- ggplot(allpairs) + 
	geom_point(aes(x=pos,y=anova_pval), size=0.5, colour="black") + 
	# geom_point(aes(x=pos,y=g1), size=0.5, colour="red") +
	# geom_point(aes(x=pos,y=g2), size=0.5, colour="blue") +
	# geom_point(aes(x=pos,y=g3), size=0.5, colour="green") +
	# geom_point(aes(x=pos,y=g4), size=0.5, colour="purple") +
	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Model with Local Ancestry') +
	theme(axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position="none") +
	geom_hline(yintercept=bonthresh)
print(g)

g <- ggplot(allpairs) + 
	geom_point(aes(x=pos,y=afr_pval), size=0.5, colour="black") + 
	# geom_point(aes(x=pos,y=g1), size=0.5, colour="red") +
	# geom_point(aes(x=pos,y=g2), size=0.5, colour="blue") +
	# geom_point(aes(x=pos,y=g3), size=0.5, colour="green") +
	# geom_point(aes(x=pos,y=g4), size=0.5, colour="purple") +
	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of AFR Local Ancestry') +
	theme(axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position="none") +
	geom_hline(yintercept=bonthresh)
print(g)

g <- ggplot(allpairs) + 
	geom_point(aes(x=pos,y=local_pval), size=0.5, colour="black") + 
	# geom_point(aes(x=pos,y=g1), size=0.5, colour="red") +
	# geom_point(aes(x=pos,y=g2), size=0.5, colour="blue") +
	# geom_point(aes(x=pos,y=g3), size=0.5, colour="green") +
	# geom_point(aes(x=pos,y=g4), size=0.5, colour="purple") +
	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Genotype with Local Ancestry') +
	theme(axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position="none") +
	geom_hline(yintercept=bonthresh)
print(g)

g <- ggplot(allpairs) + 
	geom_point(aes(x=pos,y=global_pval), size=0.5, colour="black") + 
	# geom_point(aes(x=pos,y=g1), size=0.5, colour="red") +
	# geom_point(aes(x=pos,y=g2), size=0.5, colour="blue") +
	# geom_point(aes(x=pos,y=g3), size=0.5, colour="green") +
	# geom_point(aes(x=pos,y=g4), size=0.5, colour="purple") +
	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Genotype with Global Ancestry') +
	theme(axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position="none") +
	geom_hline(yintercept=bonthresh)
print(g)

dev.off()
