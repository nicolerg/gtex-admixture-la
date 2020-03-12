# make a manhattan plot of nominal p-values for all gene-variant pairs from the local regression on chr21

library(data.table)
library(ggplot2)

# outdir <- "/users/nicolerg/gtex-admix/fixed"
# localp <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/wholeblood-chr21-local-sigdifflabels.tsv"
# gtexp <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/wholeblood-chr21-gtex-sigdifflabels.tsv"
# chr_pairs <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/wholeblood-local-allpairs-chr21.tsv.gz"

outdir <- "/users/nicolerg/gtex-admix/lung"
chr_pairs <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/lung-local-allpairs-chr1-filt.tsv.gz"
localp <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/lung-chr1-local-sigdifflabels.tsv"
gtexp <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/lung-chr1-gtex-sigdifflabels.tsv"

allpairs <- fread(input=sprintf("zcat %s", chr_pairs), sep='\t', header=TRUE, stringsAsFactors=FALSE)

bonthresh <- -log10(0.05/nrow(allpairs))
print(bonthresh)

allpairs$pos <- as.numeric(unlist(lapply(strsplit(allpairs$variant_id,split='_'), '[', 2)))
allpairs$afr_pval_nominal <- -log10(allpairs$afr_pval_nominal)
allpairs$asn_pval_nominal <- -log10(allpairs$asn_pval_nominal)
allpairs$pval_nominal <- -log10(allpairs$pval_nominal)

allpairs <- data.frame(allpairs, stringsAsFactors=FALSE)
allpairs <- allpairs[complete.cases(allpairs),]

allpairs_asn <- allpairs[order(allpairs$asn_pval_nominal, decreasing=TRUE),]
asn_sig <- unique(allpairs_asn[1:10000,1])[1:4]
asn_sig

system("echo")

allpairs_afr <- allpairs[order(allpairs$afr_pval_nominal, decreasing=TRUE),]
afr_sig <- unique(allpairs_afr[1:10000,1])[1:4]
afr_sig

afr_sig <- c("ENSG00000223745.7",
	"ENSG00000281741.2",
	"ENSG00000134202.10",
	"ENSG00000272030.1")

system("echo")

allpairs_eur <- allpairs[order(allpairs$pval_nominal, decreasing=TRUE),]
eur_sig <- unique(allpairs_eur[1:1000,1])[1:4]
eur_sig

#********************************************************************************************************
# # COLOR MOST SIGNIFICANT GENES IN EACH POPULATION

# # color significant genes in ASN

# allpairs_asn$g1 <- NA
# allpairs_asn$g2 <- NA
# allpairs_asn$g3 <- NA
# allpairs_asn$g4 <- NA

# allpairs_asn$g1[grep(asn_sig[1], allpairs_asn$gene_id, fixed=TRUE)] <- allpairs_asn$asn_pval_nominal[grep(asn_sig[1], allpairs_asn$gene_id, fixed=TRUE)]
# allpairs_asn$g2[grep(asn_sig[2], allpairs_asn$gene_id, fixed=TRUE)] <- allpairs_asn$asn_pval_nominal[grep(asn_sig[2], allpairs_asn$gene_id, fixed=TRUE)]
# allpairs_asn$g3[grep(asn_sig[3], allpairs_asn$gene_id)] <- allpairs_asn$asn_pval_nominal[grep(asn_sig[3], allpairs_asn$gene_id)]
# allpairs_asn$g4[grep(asn_sig[4], allpairs_asn$gene_id)] <- allpairs_asn$asn_pval_nominal[grep(asn_sig[4], allpairs_asn$gene_id)]

# # allpairs_asn$kcne1 <- NA
# # allpairs_asn$kcne1[grep("ENSG00000180509", allpairs_asn$gene_id)] <- allpairs_asn$asn_pval_nominal[grep("ENSG00000180509", allpairs_asn$gene_id)]
# # allpairs_asn$ybey <- NA
# # allpairs_asn$ybey[grep("ENSG00000182362", allpairs_asn$gene_id)] <- allpairs_asn$asn_pval_nominal[grep("ENSG00000182362", allpairs_asn$gene_id)]
# # allpairs_asn$sh3bgr <- NA
# # allpairs_asn$sh3bgr[grep("ENSG00000185437", allpairs_asn$gene_id)] <- allpairs_asn$asn_pval_nominal[grep("ENSG00000185437", allpairs_asn$gene_id)]
# # allpairs_asn$cct8 <- NA
# # allpairs_asn$cct8[grep("ENSG00000156261.12", allpairs_asn$gene_id)] <- allpairs_asn$asn_pval_nominal[grep("ENSG00000156261.12", allpairs_asn$gene_id)]

# jpeg(paste(outdir,'/manhattan-wholeblood-asn-siggenes.jpg',sep=''), width=1000, height=600, res=150)
# ggplot(allpairs_asn) + 
# 	geom_point(aes(x=pos,y=asn_pval_nominal), size=0.5, colour="black") + 
# 	geom_point(aes(x=pos,y=g1), size=0.5, colour="red") +
# 	geom_point(aes(x=pos,y=g2), size=0.5, colour="blue") +
# 	geom_point(aes(x=pos,y=g3), size=0.5, colour="green") +
# 	geom_point(aes(x=pos,y=g4), size=0.5, colour="purple") +
# 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Asian Local Ancestry') +
# 	theme(axis.ticks.x = element_blank(),
#       axis.text.x = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position="none")
# dev.off()

# # color significant genes in AFR

allpairs_afr$g1 <- NA
allpairs_afr$g2 <- NA
allpairs_afr$g3 <- NA
allpairs_afr$g4 <- NA

allpairs_afr$g1[grep(afr_sig[1], allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep(afr_sig[1], allpairs_afr$gene_id)]
allpairs_afr$g2[grep(afr_sig[2], allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep(afr_sig[2], allpairs_afr$gene_id)]
allpairs_afr$g3[grep(afr_sig[3], allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep(afr_sig[3], allpairs_afr$gene_id)]
allpairs_afr$g4[grep(afr_sig[4], allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep(afr_sig[4], allpairs_afr$gene_id)]

allpairs_afr$g1[allpairs_afr$g1 < bonthresh] <- NA
allpairs_afr$g2[allpairs_afr$g2 < bonthresh] <- NA
allpairs_afr$g3[allpairs_afr$g3 < bonthresh] <- NA
allpairs_afr$g4[allpairs_afr$g4 < bonthresh] <- NA

# allpairs_afr <- allpairs[order(allpairs$afr_pval_nominal, decreasing=TRUE),]
# allpairs_afr$trpm2 <- NA
# allpairs_afr$trpm2[grep("ENSG00000142185.16", allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep("ENSG00000142185.16", allpairs_afr$gene_id)]
# allpairs_afr$spatc1l <- NA
# allpairs_afr$spatc1l[grep("ENSG00000160284.14", allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep("ENSG00000160284.14", allpairs_afr$gene_id)]
# allpairs_afr$C21orf33 <- NA
# allpairs_afr$C21orf33[grep("ENSG00000160221.16", allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep("ENSG00000160221.16", allpairs_afr$gene_id)]
# allpairs_afr$pwp2 <- NA
# allpairs_afr$pwp2[grep("ENSG00000241945.7", allpairs_afr$gene_id)] <- allpairs_afr$afr_pval_nominal[grep("ENSG00000241945.7", allpairs_afr$gene_id)]

jpeg(paste(outdir,'/manhattan-wholeblood-afr-siggenes-adj.jpg',sep=''), width=1000, height=600, res=150)
ggplot(allpairs_afr) + 
	geom_point(aes(x=pos,y=afr_pval_nominal), size=0.5, colour="black") + 
	geom_point(aes(x=pos,y=g1), size=0.5, colour="red") +
	geom_point(aes(x=pos,y=g2), size=0.5, colour="blue") +
	geom_point(aes(x=pos,y=g3), size=0.5, colour="green") +
	geom_point(aes(x=pos,y=g4), size=0.5, colour="purple") +
	geom_hline(yintercept=bonthresh) +
	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of African Local Ancestry') +
	theme(axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position="none")
dev.off()

# # color significant genes in geno

# allpairs_eur$g1 <- NA
# allpairs_eur$g2 <- NA
# allpairs_eur$g3 <- NA
# allpairs_eur$g4 <- NA

# allpairs_eur$g1[grep(eur_sig[1], allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep(eur_sig[1], allpairs_eur$gene_id)]
# allpairs_eur$g2[grep(eur_sig[2], allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep(eur_sig[2], allpairs_eur$gene_id)]
# allpairs_eur$g3[grep(eur_sig[3], allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep(eur_sig[3], allpairs_eur$gene_id)]
# allpairs_eur$g4[grep(eur_sig[4], allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep(eur_sig[4], allpairs_eur$gene_id)]

# # allpairs_eur <- allpairs[order(allpairs$pval_nominal, decreasing=TRUE),]
# # allpairs_eur$pwp2 <- NA
# # allpairs_eur$pwp2[grep("ENSG00000241945", allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep("ENSG00000241945", allpairs_eur$gene_id)]
# # allpairs_eur$c21orf33 <- NA
# # allpairs_eur$c21orf33[grep("ENSG00000160221", allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep("ENSG00000160221", allpairs_eur$gene_id)]
# # allpairs_eur$umodl <- NA
# # allpairs_eur$umodl[grep("ENSG00000184385", allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep("ENSG00000184385", allpairs_eur$gene_id)]
# # allpairs_eur$fam3b <- NA
# # allpairs_eur$fam3b[grep("ENSG00000183844.16", allpairs_eur$gene_id)] <- allpairs_eur$pval_nominal[grep("ENSG00000183844.16", allpairs_eur$gene_id)]

# jpeg(paste(outdir,'/manhattan-geno-siggenes.jpg',sep=''), width=1000, height=600, res=150)
# ggplot(allpairs_eur) + 
# 	geom_point(aes(x=pos,y=pval_nominal), size=0.5, colour="black") + 
# 	geom_point(aes(x=pos,y=g1), size=0.5, colour="red") +
# 	geom_point(aes(x=pos,y=g2), size=0.5, colour="blue") +
# 	geom_point(aes(x=pos,y=g3), size=0.5, colour="green") +
# 	geom_point(aes(x=pos,y=g4), size=0.5, colour="purple") +
# 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Genotype Dosage') +
# 	theme(axis.ticks.x = element_blank(),
#       axis.text.x = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position="none")
# dev.off()

# # ********************************************************************************************************
# # COLOR VARIANT/GENE PAIRS WHOSE GENOTYPE P-VAL WAS SIGNIFICANTLY LOWER WHEN INCLUDING LOCAL COVARIATES

# sigpairs <- fread(localp, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# allpairs <- data.table(allpairs)

# setkey(allpairs, gene_id, variant_id)
# setkey(sigpairs, gene_id, variant_id)
# allpairs[sigpairs, mark := sigdiff]
# colnames(allpairs)
# print(head(allpairs))

# allpairs <- data.frame(allpairs, stringsAsFactors=FALSE)

# allpairs_asn <- allpairs[order(allpairs$asn_pval_nominal, decreasing=TRUE),]
# allpairs_afr <- allpairs[order(allpairs$afr_pval_nominal, decreasing=TRUE),]
# allpairs_eur <- allpairs[order(allpairs$pval_nominal, decreasing=TRUE),]

# allpairs_asn$s <- NA
# allpairs_asn$s[allpairs_asn$mark==1] <- allpairs_asn$asn_pval_nominal[allpairs_asn$mark==1]
# allpairs_afr$s <- NA
# allpairs_afr$s[allpairs_afr$mark==1] <- allpairs_afr$afr_pval_nominal[allpairs_afr$mark==1]
# allpairs_eur$s <- NA
# allpairs_eur$s[allpairs_eur$mark==1] <- allpairs_eur$pval_nominal[allpairs_eur$mark==1]

# jpeg(paste(outdir,'/manhattan-geno-sigdiff.jpg',sep=''), width=1000, height=600, res=150)
# ggplot(allpairs_eur) + 
# 	geom_point(aes(x=pos,y=pval_nominal), colour="black", size=0.5, alpha=0.4) + 
# 	geom_point(aes(x=pos,y=s), colour="red", size=0.5, alpha=1) + 
# 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Genotype Dosage') +
# 	theme(axis.ticks.x = element_blank(),
#       axis.text.x = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position="none") +
# 	scale_colour_manual(values=c("black","red"))
# dev.off()

# # jpeg(paste(outdir,'/manhattan-geno-sigdiff-cut.jpg',sep=''), width=1000, height=600, res=150)
# # ggplot(allpairs_eur) + 
# # 	geom_point(aes(x=pos,y=pval_nominal), colour="black", size=0.5, alpha=0.4) + 
# # 	geom_point(aes(x=pos,y=s), colour="red", size=0.5, alpha=1) + 
# # 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Genotype Dosage') +
# # 	theme(axis.ticks.x = element_blank(),
# #       axis.text.x = element_blank(),
# #       plot.title = element_text(hjust = 0.5),
# #       legend.position="none") +
# # 	scale_colour_manual(values=c("black","red")) + ylim(0,40)
# # dev.off()

# jpeg(paste(outdir,'/manhattan-asn-sigdiff.jpg',sep=''), width=1000, height=600, res=150)
# ggplot(allpairs_asn) + 
# 	geom_point(aes(x=pos,y=asn_pval_nominal), colour="black", size=0.5, alpha=0.4) + 
# 	geom_point(aes(x=pos,y=s), colour="red", size=0.5, alpha=1) + 
# 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Asian Local Ancestry') +
# 	theme(axis.ticks.x = element_blank(),
#       axis.text.x = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position="none") +
# 	scale_colour_manual(values=c("black","red"))
# dev.off()

# jpeg(paste(outdir,'/manhattan-afr-sigdiff.jpg',sep=''), width=1000, height=600, res=150)
# ggplot(allpairs_afr) + 
# 	geom_point(aes(x=pos,y=afr_pval_nominal), colour="black", size=0.5, alpha=0.4) + 
# 	geom_point(aes(x=pos,y=s), colour="red", size=0.5, alpha=1) + 
# 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of African Local Ancestry') +
# 	theme(axis.ticks.x = element_blank(),
#       axis.text.x = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position="none") +
# 	scale_colour_manual(values=c("black","red"))
# dev.off()

# #********************************************************************************************************
# # MAKE A PLOT OF THE CORRESPONDING GTEX V8 P-VALUES

# gtexpairs <- fread(gtexp, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# gtexpairs$pos <- unlist(lapply(strsplit(gtexpairs$variant_id,split='_'), '[', 2))
# gtexpairs$pval_nominal <- -log10(gtexpairs$pval_nominal)

# gtexpairs <- gtexpairs[complete.cases(gtexpairs),]

# gtexpairs$s <- NA
# gtexpairs$s[gtexpairs$sigdiff==1] <- gtexpairs$pval_nominal[gtexpairs$sigdiff==1]

# jpeg(paste(outdir,'/manhattan-gtex-sigdiff.jpg',sep=''), width=1000, height=600, res=150)
# ggplot(gtexpairs) + 
# 	geom_point(aes(x=pos,y=pval_nominal), colour="black", size=0.5, alpha=0.4) + 
# 	geom_point(aes(x=pos,y=s), colour="red", size=0.5, alpha=1) + 
# 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of GTEx Genotype Dosage') +
# 	theme(axis.ticks.x = element_blank(),
#       axis.text.x = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position="none") +
# 	scale_colour_manual(values=c("black","red"))
# dev.off()

# #********************************************************************************************************
# # MIRROR MANHATTAN: GTEX AND LAVA PVAL

# # merge allpairs_eur and gtexpairs
# gtexpairs$s <- -gtexpairs$s
# gtexpairs$pval_nominal <- -gtexpairs$pval_nominal

# colnames(gtexpairs)
# colnames(allpairs_eur)

# merged <- merge(allpairs_eur,gtexpairs, by.x=c("gene_id","variant_id","pos","mark"),by.y=c("gene_id","variant_id","pos","sigdiff"))
# head(merged)

# jpeg(paste(outdir,'/manhattan-geno-mirror.jpg',sep=''), width=1000, height=600, res=150)
# ggplot(merged) + 
# 	geom_point(aes(x=pos,y=pval_nominal.x), colour="black", size=0.5, alpha=0.4) + 
# 	geom_point(aes(x=pos,y=pval_nominal.y), colour="black", size=0.5, alpha=0.4) + 
# 	geom_point(aes(x=pos,y=s.x), colour="red", size=0.5, alpha=1) + 
# 	geom_point(aes(x=pos,y=s.y), colour="orange", size=0.5, alpha=1) + 
# 	labs(y='Nominal P-Value (-log10)', x="Position on Chromosome", title='Significance of Genotype Dosage') +
# 	theme(axis.ticks.x = element_blank(),
#       axis.text.x = element_blank(),
#       plot.title = element_text(hjust = 0.5),
#       legend.position="none") 
# dev.off()
