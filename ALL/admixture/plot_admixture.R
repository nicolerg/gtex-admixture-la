library(data.table)
library(ggplot2)
library(ggcorrplot)

### same reference panels as RFMix
admix <- fread('~/Documents/Surya/gtex-admix/admixture/merged_hg38_filt_MAF10_nomissing_pruned_original-reference-panels.3.Q.subj',header=F)
colnames(admix) <- c('SUBJID','SUB_POP','SUPER_POP','AFR','EUR','EAS')
pcs <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',sep='\t',header=T)
rfmix <- fread('~/Documents/Surya/gtex-admix/metadata/gtex_ai_global_all.txt',sep='\t',header=T)

m <- merge(admix, rfmix, by.x='SUBJID',by.y='ID', suffixes=c('_admixture','_rfmix'))

out <- c('GTEX-132QS', 'GTEX-1LGRB', 'GTEX-16AAH')
m[, color := ifelse(SUBJID %in% out, 'red','black')]

pdf('~/Desktop/admixture-v-rfmix.pdf',width=6,height=5)
ggplot(m, aes(x=AFR_admixture, y=AFR_rfmix, colour=color))+
  geom_abline(linetype='dashed') +
  geom_point()+
  scale_color_identity()

ggplot(m, aes(x=EUR_admixture, y=EUR_rfmix, colour=color))+
  geom_abline(linetype='dashed') +
  geom_point()+
  scale_color_identity()

ggplot(m, aes(x=EAS, y=ASN, colour=color))+
  geom_abline(linetype='dashed') +
  geom_point()+
  scale_color_identity()
dev.off()

gtex_meta <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)

gtex_meta <- gtex_meta[,.(SUBJID, RACE)]

pcs[,SUBJID:=sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
pcs[,c('FID','IID'):=NULL]
m <- merge(gtex_meta, pcs, by='SUBJID')
colnames(rfmix) <- paste0(colnames(rfmix),'_rfmix')
m <- merge(m, rfmix, by.x='SUBJID',by.y='ID_rfmix')
admix_gtex <- admix[grepl('GTEX',SUBJID)]
#admix_gtex[,REF_POP := NULL]
m <- merge(m, admix_gtex, by='SUBJID')

m <- m[order(EUR, AFR)]
melt <- melt(m, measure.vars=colnames(admix)[4:ncol(admix)])
pdf('~/Desktop/admixture-3-superpop-mixed.pdf',width=10,height=4)
ggplot(melt, aes(x=SUBJID,fill=variable)) +
  geom_bar(stat='identity',aes(y=value)) +
  scale_fill_manual(values=c(EUR='blue',AFR='orange',ASN='green3',EAS='green3',SAS='purple'),name='ADMIXTURE') +
  scale_x_discrete(limits=m[,SUBJID])+
  theme(axis.text.x=element_blank()) +
  labs(y='Ancestry') +
  geom_point(aes(y=-0.1,colour=factor(RACE))) +
  scale_colour_manual(values=c('1'="green3", '2'="orange", '3'="blue", '4'='#FF3D3D', '99'='gray'),
                      labels=c('1'="ASN", '2'="AFR", '3'="EUR", '4'='Native Am.', '99'='Unknown'),
                      breaks=c('3','2','1','99','4'),name='Self-Reported') 
dev.off()

m[,SUB_POP := NULL]
m[,SUPER_POP := NULL]
pdf('~/Desktop/admixture-3-superpop-mixed-corr.pdf',width=8,height=8)
ggcorrplot(cor(m[,c(3:7,23:ncol(m)),with=F]),lab=T)
dev.off()

##############################################################################################################################

# look at correlation between different versions of 3-way ADMIXTURE
admix1 <- fread('~/Documents/Surya/gtex-admix/admixture/merged_hg38_filt_MAF10_nomissing_pruned_mixed-reference-panels.3.Q.subj',header=F)
colnames(admix1) <- c('SUBJID','SUB_POP','SUPER_POP','EUR_mixed','EAS_mixed','AFR_mixed')
admix2 <- fread('~/Documents/Surya/gtex-admix/admixture/merged_hg38_filt_MAF10_nomissing_pruned_original-reference-panels.3.Q.subj',header=F)
colnames(admix2) <- c('SUBJID','SUB_POP','SUPER_POP','AFR_all-cont','EUR_CEU','EAS_CHB-JPT')
pcs <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',sep='\t',header=T)

pcs[,SUBJID:=sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
pcs[,c('FID','IID'):=NULL]

m <- merge(admix1, admix2, by=c('SUBJID','SUB_POP','SUPER_POP'))
m <- merge(m, pcs, by='SUBJID')

m[,c('SUB_POP','SUBJID','SUPER_POP') := NULL]
pdf('~/Desktop/admixture-3-superpop-comparison.pdf',width=8,height=8)
ggcorrplot(cor(m[,c(1:11),with=F]),lab=T)
ggcorrplot(cor(m[,c(1:11),with=F])[c(1:3,7:9),4:9],lab=T)
ggcorrplot(cor(m[,c(1:11),with=F])[1:3,4:6],lab=T)
dev.off()

##############################################################################################################################

# look at correlation of 2 runs of RFMix

rfmix_old <- fread('~/Documents/Surya/gtex-admix/metadata/gtex_ai_global_all.txt',sep='\t',header=T)
rfmix_new <- fread('~/Documents/Surya/gtex-admix/metadata/new_gtex_ai_global.txt',sep='\t',header=T)

m <- merge(rfmix_old, rfmix_new, by='ID', suffixes=c('_old','_new'))
pdf('~/Desktop/old-v-new-rfmix.pdf',width=6,height=5)
ggplot(m, aes(x=AFR_old, y=AFR_new))+
  geom_abline(linetype='dashed') +
  geom_point()

ggplot(m, aes(x=EUR_old, y=EUR_new))+
  geom_abline(linetype='dashed') +
  geom_point()

ggplot(m, aes(x=ASN_old, y=ASN_new))+
  geom_abline(linetype='dashed') +
  geom_point()
dev.off()



# look at all individuals 

gtex_meta <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
tg <- fread('~/Documents/Surya/gtex-admix/metadata/1kg-subject-popcodes.tsv',sep='\t',header=T)

gtex_meta <- gtex_meta[,.(SUBJID, RACE)]
tg <- tg[,.(Subject, Superpop)]
colnames(tg) <- c('SUBJID','RACE')
gtex_meta[,RACE:= as.character(RACE)]
gtex_meta[RACE=='1',RACE:='gtex_ASN']
gtex_meta[RACE=='2',RACE:='gtex_AFR']
gtex_meta[RACE=='3',RACE:='gtex_EUR']
gtex_meta[RACE=='4',RACE:='gtex_NAM']
gtex_meta[RACE=='99',RACE:='gtex_UNK']

known_pops <- rbind(gtex_meta, tg)
admix <- fread('~/Documents/Surya/gtex-admix/admixture/merged-hg38-thinned-filtered-rfmix-pruned.3.Q.subj',sep='\t',header=F)
colnames(admix) <- c('SUBJID','SUB_POP','SUPER_POP','AFR','EUR','ASN')
m <- merge(known_pops, admix, by='SUBJID')
#m[,REF_POP := NULL]

m <- m[order(EUR, AFR, ASN)]
melt <- melt(m, measure.vars=c('EUR','AFR','ASN'))
ggplot(melt, aes(x=SUBJID,fill=variable)) +
  geom_bar(stat='identity',aes(y=value)) +
  scale_fill_manual(values=c(EUR='blue',AFR='orange',ASN='green3'),name='ADMIXTURE') +
  scale_x_discrete(limits=m[,SUBJID])+
  theme(axis.text.x=element_blank()) +
  labs(y='Ancestry') +
  geom_point(data=melt[!grepl('GTEX',SUBJID)],aes(y=-0.1,colour=factor(RACE))) +
  geom_point(data=melt[grepl('GTEX',SUBJID)],aes(y=-0.05,colour=factor(RACE))) +
  scale_colour_manual(values=c(EUR='blue',gtex_EUR='blue',AMR='darkred',AFR='orange',gtex_AFR='orange',gtex_NAM='red',gtex_ASN='green3',
                               SAS='purple',EAS='green3',gtex_UNK='gray'),name='Reported (GTEx or 1000 Genomes)')





##############################################################################################################################

admix <- fread('~/Documents/Surya/gtex-admix/admixture/merged-hg38-thinned-filtered-pruned.4.Q.subj-noAMR',sep='\t',header=F)
colnames(admix) <- c('SUBJID','REF_POP','EUR','EAS','AFR','SAS')
pcs <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',sep='\t',header=T)
rfmix <- fread('~/Documents/Surya/gtex-admix/metadata/gtex_ai_global_all.txt',sep='\t',header=T)
gtex_meta <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)

gtex_meta <- gtex_meta[,.(SUBJID, RACE)]

pcs[,SUBJID:=sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
pcs[,c('FID','IID'):=NULL]
m <- merge(gtex_meta, pcs, by='SUBJID')
colnames(rfmix) <- paste0(colnames(rfmix),'_rfmix')
m <- merge(m, rfmix, by.x='SUBJID',by.y='ID_rfmix')
admix_gtex <- admix[grepl('GTEX',SUBJID)]
admix_gtex[,REF_POP := NULL]
m <- merge(m, admix_gtex, by='SUBJID')

m <- m[order(EUR, AFR, SAS, EAS)]
melt <- melt(m, measure.vars=c('EUR','EAS','AFR','SAS'))
ggplot(melt, aes(x=SUBJID,fill=variable)) +
  geom_bar(stat='identity',aes(y=value)) +
  scale_fill_manual(values=c(EUR='blue',AFR='orange',SAS='purple',EAS='green3'),name='ADMIXTURE') +
  scale_x_discrete(limits=m[,SUBJID])+
  theme(axis.text.x=element_blank()) +
  labs(y='Ancestry') +
  geom_point(aes(y=-0.1,colour=factor(RACE))) +
  scale_colour_manual(values=c('1'="green3", '2'="orange", '3'="blue", '4'='#FF3D3D', '99'='gray'),
                      labels=c('1'="ASN", '2'="AFR", '3'="EUR", '4'='Native Am.', '99'='Unknown'),
                      breaks=c('3','2','1','99','4'),name='Self-Reported') 

ggcorrplot(cor(m[,c(3:7,23:ncol(m)),with=F]))

######################
# check that reference populations are correctly assigned
admix <- fread('~/Documents/Surya/gtex-admix/admixture/merged-hg38-thinned-filtered-pruned.4.Q.subj-noAMR',sep='\t',header=F)
colnames(admix) <- c('SUBJID','REF_POP','EUR','EAS','AFR','SAS')
m <- merge(admix, tg, by.x='SUBJID',by.y='Subject')
m[,REF_POP]==m[,Superpop]

######################
# look at all individuals 

gtex_meta <- fread('~/Documents/Surya/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
tg <- fread('~/Documents/Surya/gtex-admix/metadata/1kg-subject-popcodes.tsv',sep='\t',header=T)

gtex_meta <- gtex_meta[,.(SUBJID, RACE)]
tg <- tg[,.(Subject, Superpop)]
colnames(tg) <- c('SUBJID','RACE')
gtex_meta[,RACE:= as.character(RACE)]
gtex_meta[RACE=='1',RACE:='gtex_ASN']
gtex_meta[RACE=='2',RACE:='gtex_AFR']
gtex_meta[RACE=='3',RACE:='gtex_EUR']
gtex_meta[RACE=='4',RACE:='gtex_NAM']
gtex_meta[RACE=='99',RACE:='gtex_UNK']

known_pops <- rbind(gtex_meta, tg)
admix <- fread('~/Documents/Surya/gtex-admix/admixture/merged-hg38-thinned-filtered-pruned.4.Q.subj-noAMR',sep='\t',header=F)
colnames(admix) <- c('SUBJID','REF_POP','EUR','EAS','AFR','SAS')
m <- merge(known_pops, admix, by='SUBJID')
m[,REF_POP := NULL]

m <- m[order(EUR, AFR, SAS, EAS)]
melt <- melt(m, measure.vars=c('EUR','EAS','AFR','SAS'))
ggplot(melt, aes(x=SUBJID,fill=variable)) +
  geom_bar(stat='identity',aes(y=value)) +
  scale_fill_manual(values=c(EUR='blue',AFR='orange',SAS='purple',EAS='green3'),name='ADMIXTURE') +
  scale_x_discrete(limits=m[,SUBJID])+
  theme(axis.text.x=element_blank()) +
  labs(y='Ancestry') +
  geom_point(data=melt[!grepl('GTEX',SUBJID)],aes(y=-0.1,colour=factor(RACE))) +
  geom_point(data=melt[grepl('GTEX',SUBJID)],aes(y=-0.05,colour=factor(RACE))) +
  scale_colour_manual(values=c(EUR='blue',gtex_EUR='blue',AMR='darkred',AFR='orange',gtex_AFR='orange',gtex_NAM='red',gtex_ASN='green3',
                               SAS='purple',EAS='green3',gtex_UNK='gray'),name='Reported (GTEx or 1000 Genomes)')
