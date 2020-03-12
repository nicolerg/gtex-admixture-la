# make gene sets for regulatory enrichment
library(data.table)

indir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/maf-filt/'
lavaeqtl <- paste0(indir,'gtex.admix.lava.egenes.txt.gz')
globaleqtl <- paste0(indir,'gtex.admix.global.egenes.txt.gz')

lavaeqtl <- fread(input=sprintf("zcat %s",lavaeqtl), sep='\t', header=TRUE)
globaleqtl <- fread(input=sprintf("zcat %s",globaleqtl), sep='\t', header=TRUE)

lung <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/Lung.v8.egenes.txt.gz"
lung <- fread(input=sprintf("zcat %s",lung), sep='\t', header=TRUE)

lung <- lung[,list(gene_id,pval_nominal_threshold)]
lava <- merge(lavaeqtl,lung,by="gene_id")
lava <- lava[pval_nominal < pval_nominal_threshold]
write.table(lava,paste0(indir,'/lava.gtex.cutoff.eqtl.genes.tsv'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
nrow(lava)

global <- merge(globaleqtl,lung,by="gene_id")
global <- global[pval_nominal < pval_nominal_threshold]
write.table(global,paste0(indir,'/global.gtex.cutoff.eqtl.genes.tsv'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
nrow(global)

m <- merge(lava,global,by="gene_id")
nrow(m)

# geneset 1 = same genes, different lead SNP
m <- m[variant_id.x != variant_id.y]
nrow(m)
geneset <- m[,gene_id]
write.table(geneset,paste0(indir,'/gtex.diff.lead.esnp.txt'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

# geneset 3 = top 100 genes from GEUVADIS (AFR --> EUR)
afr2eur <- '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR/AFR-EUR-eqtls-merged.RData'
load(afr2eur)
allegene[,diff:=abs(slope.EUR - slope.AFR)]
allegene <- allegene[order(diff,decreasing=TRUE)]
geneset <- allegene[1:100,gene_id]
write.table(geneset,paste0(indir,'/top100.afr2eur.txt'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

# geneset 2 = top 100 genes from GEUVADIS (EUR --> AFR)
eur2afr <- '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR-into-AFR/AFR-EUR-eqtls-merged.RData'
load(eur2afr)
allegene[,diff:=abs(slope.EUR - slope.AFR)]
allegene <- allegene[order(diff,decreasing=TRUE)]
geneset <- allegene[1:100,gene_id]
write.table(geneset,paste0(indir,'/top100.eur2afr.txt'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
