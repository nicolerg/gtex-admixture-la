library(data.table)
library(ggplot2)

tissues <- c('Brain_Frontal_Cortex_BA9','Brain_Anterior_cingulate_cortex_BA24','Brain_Amygdala')

pop <- fread('/users/nicolerg/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
pop <- pop[,list(SUBJID,RACE)]
pc <- fread('/users/nicolerg/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt', header=TRUE, sep='\t')
pc <- pc[,list(FID,PC1,PC2)]
split <- function(x){
	return(paste(unlist(unname(strsplit(x,'-')))[1:2], collapse='-'))
}
pc[,FID:=sapply(FID,split)]
m <- merge(pop, pc, by.x="SUBJID", by.y="FID")

admixed <- fread('/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt',sep='\t',header=FALSE)

# outline samples in tissue 
for (tissue in tissues){

	expr <- fread(cmd=sprintf('zcat /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/%s.v8.normalized_expression.bed.gz',tissue), sep='\t', header=T, nrow=1)
	samples <- colnames(expr)[5:ncol(expr)]
	m[,shape:=ifelse(SUBJID %in% samples,1,2)]

	m <- data.table(rbind(m, data.table(SUBJID='filler',
		RACE=10,
		PC1=NA,
		PC2=NA,
		shape=1)))
	print(m)

	g <- ggplot(data=NULL, aes(x=PC1,y=PC2)) +
		geom_point(data=m,aes(fill=factor(RACE),colour=factor(RACE),shape=factor(shape),alpha=factor(shape)),size=2) + 
		scale_fill_manual(values=c('1'="#ADFF33", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='gray', '10'=NA),guide='none') +
		scale_colour_manual(values=c('1'="#ADFF33", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='gray','10'='black'),
			labels=c('1'="Asian (N=11)", '2'="African American (N=103)", '3'="European American (N=714)", '4'='Native American (N=2)', '99'='Unknown (N=8)'),
			breaks=c('3','2','1','99','4')) +
		scale_shape_manual(values=c('1'=21,'2'=16), breaks='2', limits='2', labels='117 admixed (117AX)') +
		scale_alpha_manual(values=c('1'=0.8,'2'=0.6), guide='none') +
		geom_point(data=m[shape==1],colour='black',shape=21,alpha=0.8,size=2) +
		theme_classic() +
		theme(legend.justification=c(0,0), 
			legend.position='none',
			legend.text = element_text(size=10),
			legend.title = element_blank(),
	    	legend.spacing = unit(3,unit='mm')) +
		guides(colour = guide_legend(override.aes = list(size=2,
									alpha=1,
									shape=21,
									colour=NA,
									fill=c("#0000FF","#FF9900","#ADFF33",'gray','#FF3D3D')), keywidth=5, keyheight=5, default.unit="mm"),
				shape = guide_legend(override.aes = list(size=2,
					alpha=1,
					shape=21,
					colour='black',fill='white'), keywidth=5, keyheight=5, default.unit="mm")) +
		labs(x='gPC1',y='gPC2')

	pdf(sprintf("~/gtex-admix/plots/september/geno_pc_%s.pdf",tissue),width=4,height=4)
	print(g)
	dev.off()
}
