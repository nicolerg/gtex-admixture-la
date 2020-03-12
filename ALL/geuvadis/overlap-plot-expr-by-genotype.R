library(data.table)
library(ggplot2)
library(gridExtra)

infile <- '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/overlapping/overlapping.snp.expression.tsv.gz'
dt <- fread(input=sprintf("zcat %s", infile), sep='\t', header=TRUE)

pdf("~/gtex-admix/plots/afr-eur-overlapping-eqtl.pdf",width=4,height=4)
for (snp in unique(dt[,variant])){
	sub <- dt[variant==snp]
	g <- dt[variant==snp,gene][1]
	p <- ggplot(sub,aes(x=factor(genotype),y=expr,colour=factor(pop)))+
	geom_boxplot() +
	theme_bw() + 
	scale_colour_manual(values=c(AFR='orange',EUR='blue')) +
	theme(axis.title.x=element_blank(),
		axis.text.x=element_text(size=14),
		legend.title=element_blank()) +
	labs(y="Normalized expr", title=paste0(g,':',snp)) 
	print(p)
}
dev.off()