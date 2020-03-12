library(data.table)
library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outdir <- paste0(args[2],'/')

if(grepl('AFR-into-EUR', outdir)){
	out1 <- '~/gtex-admix/plots/GEU.AFR.eqtl.0EUR.var.expr.pdf'
	out2 <- paste0(outdir,'GEU.AFR.eqtl.0EUR.tstat.tsv')
	fillm <- "blue"
	fillp <- "orange"
} else if (grepl('EUR-into-AFR', outdir)){
	out1 <- '~/gtex-admix/plots/GEU.EUR.eqtl.0AFR.var.expr.pdf'
	out2 <- paste0(outdir,'GEU.EUR.eqtl.0AFR.tstat.tsv')
	fillp <- "blue"
	fillm <- "orange"
}

tstat <- file(out2, 'w')
cat("gene_id\tvariant_id\tpoly_genotype\tmono_genotype\ttstat\n", file=tstat, append=TRUE)

dt <- fread(input=sprintf("zcat %s", infile), sep='\t', header=TRUE)

pairs <- dt[,list(gene,variant)]
pairs <- unique(pairs)

ttest <- function(mono,poly){
	# perform a t-test between AFR expression and each EUR genotype expression
	poly[,ttest := '']
	for (geno in unique(poly[,genotype])){
		max <- max(poly[genotype == geno, expr])
		pval <- unname(unlist(t.test(poly[,genotype==geno,expr],mono[,expr])['statistic']))
		string <- paste(g,var,geno,unique(mono[,genotype])[1],pval,sep='\t')
		print(string)
		cat(paste(string,'\n',sep=''), file=tstat, append=TRUE)
		poly[genotype == geno & expr == max, ttest := as.character(round(pval,3))]
	}

	poly[,genotype := as.factor(genotype)]
	mono[,genotype := as.factor(genotype)]

	# multiple groups for EUR
	p <- ggplot(data=poly, aes(x=genotype,expr)) +
		geom_boxplot(fill=fillp) +
		geom_jitter(width = 0.2, size=0.5) +
		theme_bw() +
		theme(axis.title.x=element_blank(),legend.position="none",
			axis.text.x=element_text(size=14)) +
		labs(y="Normalized expr", title=paste0(g,':',var)) +
		geom_text(aes(label=ttest,
			y = expr + 0.1),
			position = position_dodge(0.9),
			vjust = 0)
	# single genotype for AFR
	m <- ggplot(mono, aes(x=genotype,y=expr)) +
		geom_boxplot(fill=fillm) +
		geom_jitter(width = 0.2, size = 0.5) +
		theme_bw() +
		theme(axis.title=element_blank(),legend.position="none",
			axis.text.x=element_text(size=14),
			title=element_blank()) +
		coord_cartesian( ylim = c(round(min(poly[,expr]),1) , round(max(poly[,expr]),1)) )

	# print them next to each other
	return(grid.arrange(p,m,nrow=1,ncol=2,widths=c(4,2)))
}

pdf(out1,width=6,height=4)
for(i in 1:nrow(pairs)){
	g <- pairs[i,gene]
	var <- pairs[i,variant]
	sub <- dt[gene == g & variant == var]

	afr <- sub[pop == "AFR"]
	eur <- sub[pop == "EUR"]

	if(grepl('AFR-into-EUR', outdir)){
		ttest(eur,afr)
	} else if (grepl('EUR-into-AFR', outdir)){
		ttest(afr,eur)
	}
}
dev.off()

close(tstat)

