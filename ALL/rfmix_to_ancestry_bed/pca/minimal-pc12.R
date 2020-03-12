library(data.table)
library(ggplot2)
library(plyr)

infile <- '/mnt/lab_data/montgomery/nicolerg/pca-subset/gtexV8-PC12-df.tsv'
dat <- fread(infile, sep='\t', header=TRUE)

group.colors <- c(ESN="#FFCD00",
	GWD="#FFB900",
	LWK="#CC9933",
	MSL="#E1B919",
	YRI="#FFB933",
	ACB="#FF9900",
	ASW="#FF6600",
	CLM="#CC3333",
	MXL="#E10033",
	PEL="#FF0000",
	PUR="#CC3300",
	CDX="#339900",
	CHB="#ADCD00",
	CHS="#00FF00",
	JPT="#008B00",
	KHV="#00CC33",
	CEU="#0000FF",
	FIN="#00C5CD",
	GBR="#00EBFF",
	IBS="#6495ED",
	TSI="#00008B",
	BEB="#8B008B",
	GIH="#9400D3",
	ITU="#B03060",
	PJL="#E11289",
	STU="#FF00FF")

breaks <- c("ESN",
	"GWD",
	"LWK",
	"MSL",
	"YRI",
	"ACB",
	"ASW",
	"CLM",
	"MXL",
	"PEL",
	"PUR",
	"CDX",
	"CHB",
	"CHS",
	"JPT",
	"KHV",
	"CEU",
	"FIN",
	"GBR",
	"IBS",
	"TSI",
	"BEB",
	"GIH",
	"ITU",
	"PJL",
	"STU")

pdf('~/gtex-admix/plots/final_figures/gtex-TG-geno-PC12.pdf', width=8, height=6)
ggplot(data=dat) +
	geom_point(data=subset(dat, !pop=="GTEx"), aes(x=PC1, y=PC2, colour=pop, shape=dataset), alpha=1, size=1) + 
	theme_bw() + 
	theme(legend.position="right",
		plot.title=element_text(hjust = 0.5),
		legend.text=element_text(size=12),
		legend.title=element_text(size=15),
		axis.text=element_text(size=10),
		legend.key.size=unit(18,"point"),
		axis.title=element_text(size=15)) +
	scale_colour_manual(values=group.colors,name="1000 Genomes",breaks=breaks) +
	geom_point(data=subset(dat, pop=="GTEx"), aes(x=PC1, y=PC2, shape=dataset), colour="black", size=1.5) +
	scale_shape_manual(values=c(20,2),name="Consortium") +
	xlab("PC1") + ylab("PC2") 
dev.off()

pdf('~/gtex-admix/plots/final_figures/gtex-geno-PC12.pdf', width=6, height=6)
ggplot(data=dat) +
	theme_bw() + 
	geom_point(data=subset(dat, pop=="GTEx"), aes(x=PC1, y=PC2), colour="black", size=1.5, shape=2) +
	xlab("PC1") + ylab("PC2") 
dev.off()

# mock PCA plot 

mock <- dat[dataset=='GTEx']
mock <- mock[Race == 'EUR' | Race == 'AFR' | Race == 'ASN']
head(mock)

mock[PC2 < 0.02 & PC1 < 0.01, Race := 'EUR']
mock[PC2 > 0.02, Race := 'ASN']

pdf('~/gtex-admix/plots/final_figures/gtex-race-geno-PC12.pdf', width=6, height=4)
ggplot(data=mock, aes(x=PC1, y=PC2, colour=factor(Race))) +
	geom_point() + 
	theme_bw() + 
	scale_colour_manual(values=c(ASN="#ADFF33", AFR="#FF9900", EUR="#0000FF"),
		labels=c(ASN="Asian", AFR="African American", EUR="White"),
		name="Self-Reported Race",
		breaks=c('EUR','AFR','ASN')) +
	xlab("PC1") + ylab("PC2") 
dev.off()

dat <- dat[dataset=='GTEx']

# genotype PC plot with selected individuals shown 

pop <- fread('/users/nicolerg/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
pop <- pop[,list(SUBJID,RACE)]

pcs <- fread('/mnt/lab_data/montgomery/nicolerg/pca-subset/gtex-tg-intersection-pcs.csv', sep=',', header=FALSE)
head(pcs)
colnames(pcs) <- c('SUBJID','PC1','PC2','PC3','PC4')
m <- merge(dat, pcs, by=c('PC1', 'PC2'))

admixed <- fread('/mnt/lab_data/montgomery/nicolerg/rfmix-admixed/combined/gtex-admixed0.9.txt',sep='\t',header=FALSE)
m[,shape:=ifelse(SUBJID %in% admixed[,V1],1,2)]

geno_pc <- ggplot(data=NULL, aes(x=PC1,y=PC2)) +
	geom_point(data=m,aes(fill=factor(Race),colour=factor(Race),shape=factor(shape),alpha=factor(shape)),size=1) + 
	scale_fill_manual(values=c(ASN="#ADFF33", AFR="#FF9900", EUR="#0000FF", IND='#FF3D3D', OTH='gray'),
		labels=c(ASN="Asian (N=11)", AFR="African American (N=103)", EUR="White (N=714)", IND='American Indian (N=2)', OTH='Unknown (N=8)'),
		name="Self-Reported Race",guide='none') +
	scale_colour_manual(values=c(ASN="#ADFF33", AFR="#FF9900", EUR="#0000FF", IND='#FF3D3D', OTH='gray'),
		labels=c(ASN="Asian (N=11)", AFR="African American (N=103)", EUR="White (N=714)", IND='American Indian (N=2)', OTH='Unknown (N=8)'),
		name="Self-Reported Race",
		breaks=c('EUR','AFR','ASN','OTH','IND')) +
	scale_shape_manual(name="Global Ancestry",labels=c('1'='<90% from single population','2'='>=90% from single population'),values=c('1'=21,'2'=16)) +
	scale_alpha_manual(values=c('1'=1,'2'=0.7), guide='none') +
	geom_point(data=m[shape==1],colour='black',shape=21,alpha=0.8,size=1.5) +
	theme_bw() +
	#theme(legend.position=c(0.3,0.6),legend.background = element_rect(color = "black", size = 0.3, linetype = "solid")) +
	guides(colour = guide_legend(override.aes = list(size=2,alpha=1), order=1),
		shape = guide_legend(override.aes = list(size=c(2,1),alpha=1), order=0))

pdf('~/gtex-admix/plots/final_figures/gtex-race-geno-PC12-mypca.pdf', width=7, height=4)
geno_pc
dev.off()


