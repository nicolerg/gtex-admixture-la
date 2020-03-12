library(data.table)
library(ggplot2)

indir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/'
allpairs <- paste0(indir,'lung-LAVA-global-allpairs-merged.tsv.gz')

lavaeqtl <- fread(paste0(indir,'lava.eqtl.genes.tsv'), sep='\t', header=TRUE, select=c("gene_id","variant_id","qval"))
globaleqtl <- fread(paste0(indir,'global.eqtl.genes.tsv'), sep='\t', header=TRUE, select=c("gene_id","variant_id","qval"))

allp <- fread(input=sprintf("zcat %s",allpairs), sep='\t', header=TRUE, select=c('gene_id','variant_id','global_pval_nominal','global_slope','lava_pval_nominal','lava_slope'))

allp <- allp[global_pval_nominal < 0.5 & lava_pval_nominal < 0.5]

# plot pval vs pval

jpeg("~/gtex-admix/plots/gtex.lung.global.local.pval.jpg", res=200, width=1000, height=1000)
ggplot(allp, aes(x=global_pval_nominal, y=lava_pval_nominal)) +
	geom_point() +
	theme_bw() +
	labs(x = "Global nominal pvalue", y = "LAVA nominal pvalue")
dev.off()

# filter down pairs; plot effect sizes

allp <- allp[global_pval_nominal < 0.01 & lava_pval_nominal < 0.01]
eqtls <- merge(lavaeqtl, globaleqtl, by=c("gene_id","variant_id"), all=TRUE)
eqtls[,qval.y := NULL]
eqtls[,qval.x := 0]
merged <- merge(allp, eqtls, by = c("gene_id","variant_id"), all=TRUE)
merged[qval.x != 0, qval.x := 1]

jpeg("~/gtex-admix/plots/gtex.lung.global.local.slope.jpg", res=200, width=1000, height=1000)
ggplot(merged, aes(x=global_slope, y=lava_slope, colour=factor(qval.x), alpha=factor(qval.x))) +
	geom_point() +
	theme_bw() +
	scale_colour_manual(values=c('0'='red','1'='black'),guide='none') +
	scale_alpha_manul(values=c('0'=1,'1'=0.3),guide='none')
	labs(x = "Global effect size", y = "LAVA effect size")
dev.off()
