library(data.table)

outdir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/anova/filtered'

allpairsf <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/lung-local-allpairs-chr1-filt.tsv.gz'
# gene_id, variant_id, afr_pval_nominal
genof <- '/mnt/lab_data/montgomery/nicolerg/filtered/gtex.MAF01.phased.chr1.vcf.gz'
# ID = variant ID
localcovf <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr1.hg19.localcov.tsv.gz'
# SNP_ID = variant ID (multiple rows per variant)
exprf <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/whole_chr/lung.chr1.expression.bed'

geno <- fread(input=sprintf("zcat %s", genof), sep='\t', header=TRUE, stringsAsFactors=FALSE)
localcov <- fread(input=sprintf("zcat %s", localcovf), header=TRUE, sep='\t', stringsAsFactors=FALSE)
allpairs <- fread(input=sprintf("zcat %s", allpairsf), header=TRUE, sep='\t', stringsAsFactors=FALSE)
expr <- fread(exprf, header=FALSE, sep='\t', stringsAsFactors=FALSE)
subj <- paste0('subj',1:(ncol(expr)-4))
colnames(expr) <- c('chr','start','end','gene_id',subj)

# filter allpairs to include only pairs that reach significance for African local ancestry 
allpairs[,afr_pval_nominal := as.numeric(allpairs[,afr_pval_nominal])]
print('allpairs')
nrow(allpairs)
allpairs <- allpairs[ afr_pval_nominal < 3.16e-9 ] #10e-8.5
nrow(allpairs)

keepsnp <- unique(allpairs[,variant_id])
keepgene <- unique(allpairs[,gene_id])

# filter geno 
print('geno')
nrow(geno)
genofilt <- subset(geno, ID %in% keepsnp)
nrow(genofilt)
# filter localcov 
print('localcov')
nrow(localcov)
lcovfilt <- subset(localcov, SNP_ID %in% keepsnp)
nrow(lcovfilt)
# filter expression 
print('expr')
nrow(expr)
exprfilt <- subset(expr, gene_id %in% keepgene)
nrow(exprfilt)

write.table(exprfilt, file=paste0(outdir,'/chr1-expr-filt-byafrsig.tsv'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(genofilt, file=paste0(outdir,'/chr1-geno-filt-byafrsig.tsv'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(lcovfilt, file=paste0(outdir,'/chr1-localcov-filt-byafrsig.tsv'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
