library(data.table)
library(ggplot2)

het <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/het_fraction/diff-lead-pairs-by-tissue-1e04-het-frac.tsv',sep='\t',header=T)
pairs <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/het_fraction/diff-lead-pairs-by-tissue-1e-04.txt',sep='\t',header=T)

m <- merge(pairs, het, by.y=c('variant_id','gene','tissue'), by.x=c('variant_id_local','gene_id','tissue'),all=T)
m <- merge(m, het, by.y=c('variant_id','gene','tissue'), by.x=c('variant_id_global','gene_id','tissue'), suffixes=c('_local','_global'))

m <- unique(m)
head(m)

m[,het_frac_diff := het_frac_local - het_frac_global]
m[,pval_diff := -log10(pval_nominal_local) + log10(pval_nominal_global)]

head(m)

pdf('~/gtex-admix/plots/final_figures/het-frac-1e04.pdf',width=8,height=4)
ggplot(m, aes(x=het_frac_diff, y=pval_diff, colour=tissue)) +
	geom_point(alpha=0.5)

ggplot(m, aes(x=het_frac_diff, y=D, colour=tissue)) +
	geom_point(alpha=0.5)

ggplot(m, aes(x=pval_diff, y=D, colour=tissue)) +
	geom_point(alpha=0.5)


# are het fractions sig dif within each tissue?
m <- m[D < 1]
f <- m[,.(tissue, het_frac_local, het_frac_global)]
f <- melt(f, id.vars='tissue')

ggplot(f, aes(x=tissue, y=value, fill=variable))+
	geom_boxplot()

ggplot(f, aes(x=variable, y=value))+
	geom_boxplot()

dev.off()

# does heterozygosity influence differences in colocalization?



