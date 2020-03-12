library(data.table)
library(ggplot2)

pop <- fread('GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
pop <- pop[,list(SUBJID,RACE)]
pc <- fread('GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt', header=TRUE, sep='\t')
pc <- pc[,list(FID,PC1,PC2)]
split <- function(x){
  return(paste(unlist(unname(strsplit(x,'-')))[1:2], collapse='-'))
}
pc[,FID:=sapply(FID,split)]
head(pc)
head(pop)
m <- merge(pop, pc, by.x="SUBJID", by.y="FID")
head(m)
# =====
## select the most admixed individuals
# m <- m[PC1 < 0 | PC2 > 0.025]
# write.table(m,'gtex-most-admixed.tsv',sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
# sub <- m[,SUBJID]
# write.table(sub,'gtex-most-admixed-subj.txt',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)


# ===== 
m[,shape:=ifelse(RACE=='1'|RACE=='2',1,2)]

pdf('~/Documents/Surya/gtex-admix/plots/figures/GTEx_v8_genotype_pca.pdf', height=4, width=4)
ggplot(data=NULL, aes(x=PC1,y=PC2)) +
  geom_point(data=m,aes(fill=factor(RACE),colour=factor(RACE),shape=factor(shape),alpha=factor(shape)),size=0.8) + 
  scale_fill_manual(values=c('1'="#ADFF33", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='black'),
                      labels=c('1'="Asian (N=11)", '2'="African American (N=103)", '3'="White (N=714)", '4'='American Indian (N=2)', '99'='Unknown (N=8)'),
                      name="Self-Reported Race",guide='none') +
  scale_colour_manual(values=c('1'="#ADFF33", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='black'),
                    labels=c('1'="Asian (N=11)", '2'="African American (N=103)", '3'="White (N=714)", '4'='American Indian (N=2)', '99'='Unknown (N=8)'),
                    name="Self-Reported Race",
                    breaks=c('3','2','1','99','4')) +
  scale_shape_manual(name="Local Ancestry",labels=c('1'='Imputed','2'='Not imputed'),values=c('1'=21,'2'=16)) +
  scale_alpha_manual(values=c('1'=1,'2'=0.5), guide='none') +
  geom_point(data=m[RACE==1|RACE==2],colour='black',shape=21,alpha=0.8) +
  theme_bw() +
  theme(legend.position=c(0.3,0.6),legend.background = element_rect(color = "black", size = 0.3, linetype = "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2,alpha=1)),
         shape = guide_legend(override.aes = list(size=c(2,1),alpha=1)))
dev.off()
