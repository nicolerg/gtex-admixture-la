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

global[, gtex_cutoff := -log10(gtex_cutoff)]
global[, pval_nominal := -log10(pval_nominal)]
lava[, gtex_cutoff := -log10(gtex_cutoff)]
lava[, pval_nominal := -log10(pval_nominal)]

# get number of eGenes with each method
print(paste('Number of eGenes in global:',length(unique(global[,gene_id]))))
print(paste('Number of eGenes in LAVA:',length(unique(lava[,gene_id]))))

# filter down to those only present in one method 
global_uniq <- global[!(gene_id %in% lava[,gene_id])]
paste('Unique global eGenes:',length(unique(global_uniq[,gene_id])))
lava_uniq <- lava[!(gene_id %in% global[,gene_id])]
paste('Unique LAVA eGenes:',length(unique(lava_uniq[,gene_id])))

# make a plot of # eGenes versus magnitude above cutoff (total; not unique)

m <- max(
	lava[,pval_nominal] - lava[,gtex_cutoff],
	global[,pval_nominal] - global[,gtex_cutoff])

data <- data.table(mag_above=seq(2,m,0.1),
	n_global_egenes=0,
	n_lava_egenes=0)

for (cutoff in data[,mag_above]){
	lava_tmp <- lava[ (pval_nominal - gtex_cutoff) >= cutoff ]
	global_tmp <- global[ (pval_nominal - gtex_cutoff) >= cutoff ]
	data[ mag_above==cutoff, n_lava_egenes := length(unique(lava_tmp[,gene_id])) ]
	data[ mag_above==cutoff, n_global_egenes := length(unique(global_tmp[,gene_id])) ]
}

data <- melt(data,id='mag_above')

jpeg(paste0('/users/nicolerg/gtex-admix/plots/',tissue,'-egenes-per-cutoff.jpg'), width=1600, height=1000, res=200)
ggplot(data, aes(x=mag_above, y=value, colour=factor(variable))) + 
	geom_line() +
	geom_point() +
	labs(y='N eGenes', x='Orders of magnitude above p-value cutoff', title=tissue) +
	scale_colour_manual(values=c(n_global_egenes='blue',n_lava_egenes='red'),
		labels=c(n_global_egenes='Global',n_lava_egenes='LAVA'),name='Method')
dev.off()

# now UNIQUE

m <- max(
	lava[,pval_nominal] - lava[,gtex_cutoff],
	global[,pval_nominal] - global[,gtex_cutoff])

data2 <- data.table(mag_above=seq(1,m,0.1),
	n_global_egenes=0,
	n_lava_egenes=0)

for (cutoff in data2[,mag_above]){
	lava_less_stringent <- lava[ (pval_nominal - gtex_cutoff) >= cutoff - 1 ]
	global_less_stringent <- global[ (pval_nominal - gtex_cutoff) >= cutoff - 1 ]

	lava_tmp <- lava[ (pval_nominal - gtex_cutoff) >= cutoff ]
	global_tmp <- global[ (pval_nominal - gtex_cutoff) >= cutoff ]
	
	g_uniq <- global_tmp[!(gene_id %in% lava_less_stringent[,gene_id])]
	l_uniq <- lava_tmp[!(gene_id %in% global_less_stringent[,gene_id])]
	
	data2[ mag_above==cutoff, n_lava_egenes := length(unique(l_uniq[,gene_id])) ]
	data2[ mag_above==cutoff, n_global_egenes := length(unique(g_uniq[,gene_id])) ]
}

data2 <- melt(data2,id='mag_above')

jpeg(paste0('/users/nicolerg/gtex-admix/plots/',tissue,'-uniq-egenes-per-cutoff.jpg'), width=1600, height=1000, res=200)
ggplot(data2, aes(x=mag_above, y=value, colour=factor(variable))) + 
	geom_line() +
	geom_point() +
	labs(y='N eGenes', x='Magnitude above p-value cutoff (-log10)', title=paste(tissue,'(unique)')) +
	scale_colour_manual(values=c(n_global_egenes='blue',n_lava_egenes='red'),
		labels=c(n_global_egenes='Global',n_lava_egenes='LAVA'),name='Method')
dev.off()

# write out files of unique-to-method eSNPs to annotate

# limit it to eGenes where pval in one method is >= 1 order of magnitude larger than in other method 

lava_high <- lava[ (pval_nominal - gtex_cutoff) >= 1 ]
global_high <- global[ (pval_nominal - gtex_cutoff) >= 1 ]

lava_uniq <- lava_high[ !(gene_id %in% global[,gene_id]) ]
global_uniq <- global_high[ !(gene_id %in% lava[,gene_id]) ]

lava_uniq[,method:='LAVA']
global_uniq[,method:='global']
master <- data.table(rbind(lava_uniq,global_uniq))


paste(tissue,'lava',length(unique(master[method=='LAVA',gene_id])))
paste(tissue,'global',length(unique(master[method=='global',gene_id])))

# write out these files for annotation
write.table(master, paste0(outdir,'/',tissue,'-uniq-eqtls-ties.txt'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
