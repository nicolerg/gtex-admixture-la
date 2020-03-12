library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outdir <- args[2]

eur2afr <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR-into-AFR/tstat-per-genotype.tsv"
afr2eur <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR/tstat-per-genotype.tsv"

eur2afr <- fread(eur2afr, sep='\t', header=TRUE)
afr2eur <- fread(afr2eur, sep='\t', header=TRUE)

outfile <- "~/gtex-admix/plots/figures/geu-tstat.pdf"

eur2afr[,pop:="EUR"]
afr2eur[,pop:="AFR"]

all <- rbind(eur2afr,afr2eur)

pdf(outfile,width=4,height=4)
ggplot(all, aes(x=factor(label),y=tstat,colour=factor(pop))) +
	geom_boxplot() +
	theme_bw() +
	scale_x_discrete(labels=c(altalt="Homozygous\ndifferent",refref="Homozygous\nsame",het='Heterozygous')) +
	scale_colour_manual(name="",values=c(EUR="blue",AFR="orange"),labels=c(EUR="EUR into AFR",AFR="AFR into EUR")) +
	labs(y="T-test statistic") +
	theme(axis.title.x=element_blank(),
		axis.text.x=element_text(size=10,colour="black"),
		legend.title=element_blank(),
		legend.position=c(0.81,0.9),
		legend.background = element_rect(color = "black", size = 0.3, linetype = "solid")) +
	ylim(c(-20,20))
dev.off()
