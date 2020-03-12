library(data.table)
library(ggcorrplot)

pcs <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',sep='\t',header=T)
global_ai <- fread('~/Documents/Surya/gtex-admix/metadata/gtex_ai_global_all.txt',sep='\t',header=T)
pcs[,SUBJID := sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
m <- merge(pcs, global_ai, by.x='SUBJID', by.y='ID')
m[,UNK := NULL]
m <- data.frame(m)
rownames(m) <- m$SUBJID
m$SUBJID <- NULL
m$FID <- NULL
m$IID <- NULL
corr <- cor(m)
ggcorrplot(corr[1:10, 21:23],lab=T)
