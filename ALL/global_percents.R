library(data.table)

gai <- fread('~/Documents/Surya/gtex-admix/metadata/gtex_ai_global_all.txt',sep='\t',header=T)
ids <- fread('~/Documents/Surya/gtex-admix/metadata/gtex-admixed0.9.txt',sep='\t',header=F)
head(gai)

gai <- gai[ID %in% ids[,V1]]
asn <- sum(gai[,ASN])/nrow(gai)
afr <- sum(gai[,AFR])/nrow(gai)
eur <- sum(gai[,EUR])/nrow(gai)
unk <- sum(gai[,UNK])/nrow(gai)

sprintf('ASN: %s',asn)
sprintf('AFR: %s',afr)
sprintf('EUR: %s',eur)
sprintf('UNK: %s',unk)