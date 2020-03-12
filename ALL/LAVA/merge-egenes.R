library(data.table)
library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly=TRUE)
global <- args[1]
lava <- args[2]
outdir <- args[3]
tissue <- args[4]

# indir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung'
# global <- paste0(indir,'/gtex.admix.global.egenes.tied.txt.gz')
# lava <- paste0(indir,'/gtex.admix.lava.egenes.tied.txt.gz')
# outdir <- indir
# tissue <- 'Lung'

global <- fread(cmd=sprintf("zcat %s",global), sep='\t', header=TRUE)
lava <- fread(cmd=sprintf("zcat %s",lava), sep='\t', header=TRUE)

print(tissue)
ngenes <- length(unique(lava[,gene_id]))
print(paste('Number of genes:',ngenes))

# filter down to those that pass the GTEx threshold
global <- global[pval_nominal < gtex_cutoff]
lava <- lava[pval_nominal < gtex_cutoff]

# write out these files for annotation
write.table(global, paste0(outdir,'/',tissue,'-global-eqtls-ties.txt'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(lava, paste0(outdir,'/',tissue,'-lava-eqtls-ties.txt'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

# get number of eGenes with each method
globale <- length(unique(global[,gene_id]))
print(paste('Number of eGenes in global:',globale))
lavae <- length(unique(lava[,gene_id]))
print(paste('Number of eGenes in LAVA:',lavae))

# get number of eSNPs from each eGene 
global_counts <- as.data.table(table(global[,gene_id]))
colnames(global_counts) <- c('gene_id','n_esnps_global')
lava_counts <- as.data.table(table(lava[,gene_id]))
colnames(lava_counts) <- c('gene_id','n_esnps_lava')
merged <- merge(global_counts, lava_counts, by='gene_id')
overlap <- length(merged[,gene_id])
print(paste('Number of overlapping eGenes:',overlap))
m <- melt(merged, id='gene_id')

# box and whisker plot
pdf(paste0('/users/nicolerg/gtex-admix/plots/',tissue,'-esnps-per-egene-boxplot.pdf'), width=6, height=4)
ggplot(m, aes(x=variable, y=value)) + 
	geom_boxplot() +
	theme(axis.title.x=element_blank()) +
	labs(y='N tied top SNPs', title='Number of tied eSNPs in matching eGenes') +
	scale_x_discrete(limits=c('n_esnps_global', 'n_esnps_lava'), 
		labels=c(n_esnps_lava='LAVA',n_esnps_global='global'))
dev.off()

# density plot?
pdf(paste0('/users/nicolerg/gtex-admix/plots/',tissue,'-esnps-per-egene-density.pdf'), width=6, height=4)
ggplot(m, aes(group=variable, x=value, fill=variable)) + 
	geom_density(alpha=0.5) +
	labs(x='N tied top SNPs', title='Number of tied eSNPs in matching eGenes') +
	scale_fill_manual(values=c(n_esnps_lava='red',n_esnps_global='blue'),
		labels=c('n_esnps_lava'='LAVA','n_esnps_global'='global'),
		name='Method') 
dev.off()

pdf(paste0('/users/nicolerg/gtex-admix/plots/',tissue,'-esnps-per-egene-density-zoom.pdf'), width=6, height=4)
ggplot(m, aes(group=variable, x=value, fill=variable)) + 
	geom_density(alpha=0.5) +
	labs(x='N tied top SNPs', title='Number of tied eSNPs in matching eGenes') +
	scale_fill_manual(values=c(n_esnps_lava='red',n_esnps_global='blue'),
		labels=c('n_esnps_lava'='LAVA','n_esnps_global'='global'),
		name='Method') +
	coord_cartesian(xlim=c(1,45),ylim=c(0,0.015))
dev.off()

# number of overlapping eGenes with different lead SNPs; write these out 
# i.e. where the length of the intersection between tied SNP sets is 0
l <- list()
i <- 1
mismatch = 0
for (gene in merged[,gene_id]){
	lava_esnps <- lava[gene_id == gene, variant_id]
	global_esnps <- global[gene_id == gene, variant_id]
	if (length(intersect(lava_esnps,global_esnps)) == 0) {
		mismatch = mismatch + 1
		ltmp <- lava[gene_id == gene]
		gtmp <- global[gene_id == gene]
		ltmp[,method := 'LAVA']
		gtmp[,method := 'global']
		m <- rbind(ltmp,gtmp)
		l[[i]] <- m
		i <- i + 1
	}
}
master <- rbindlist(l)
write.table(master, paste0(outdir, '/', tissue, '.samegene.diffSNP.master.tsv'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
print(paste('Number of matching eGenes where leads SNPs are different:',mismatch))

# compute mcnemar chi-squred test for correlated data

a <- overlap
b <- globale - a 
c <- lavae - a 
d <- ngenes - (a + b + c)

m <- matrix(c(a,b,c,d),nrow=2,ncol=2,byrow=TRUE)
t <- as.table(m)
colnames(t) <- c('eqtl','not_eqtl')
rownames(t) <- colnames(t)

mcnemar.test(t)
