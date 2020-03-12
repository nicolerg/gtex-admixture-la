library(SNPRelate)
library(data.table)
library(ggplot2)
library(plyr)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#args <- commandArgs(trailingOnly=TRUE)
outdir <- "/mnt/lab_data/montgomery/nicolerg/pca-subset"
vcf.in.1kg <- "/mnt/lab_data/montgomery/nicolerg/pca-subset/tg.master.subset.vcf"
vcf.in.gtex <- "/mnt/lab_data/montgomery/nicolerg/pca-subset/gtex.master.subset.vcf"
tg.pop.codes <- "/users/nicolerg/gtex-admix/metadata/tg_pop_labels.tsv"
tg.pop.match <- "/users/nicolerg/gtex-admix/metadata/1000genomes_match_popcodes.tsv"
tg.colors <- "/users/nicolerg/gtex-admix/metadata/1kg_pop_hex.tsv"
phenofile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

# # --------------------------------------------------------------------------------------------------------------
# # IF YOU HAVEN'T ALREADY MADE A SINGLE VCF OF THE INTERSECTION
# # read in GTEx
# system("echo 'Reading in GTEx overlap...'")
# gtex <- fread(input = sprintf("zcat %s", vcf.in.gtex), header=TRUE, sep="\t", stringsAsFactors=FALSE)
# # read in 1000 Genomes
# system("echo 'Reading in 1000 Genomes overlap...'")
# tg <- fread(input = sprintf("zcat %s", vcf.in.1kg), header=TRUE, sep="\t", stringsAsFactors=FALSE)
# # sort both of the data tables
# system("echo 'Sorting the variants...'")
# tg <- tg[order(POS)]
# gtex <- gtex[order(POS)]
# # find the intersection
# system("echo")
# system("echo 'Finding overlap between data sets'")
# print(paste("N GTEx variants:", nrow(gtex), sep=' '), quote=FALSE)
# print(paste("N 1000 Genome variants:", nrow(tg), sep=' '), quote=FALSE)
# # select only shared variants 
# shared <- intersect(gtex$POS, tg$POS)
# nshared <- length(shared)
# print(paste("N shared variants:", length(shared), sep=' '), quote=FALSE)
# gtex <- gtex[match(shared,gtex$POS),]
# tg <- tg[match(shared,tg$POS),]
# system("echo")
# # get rid of leading columns (1-9) in 1000 Genomes and cbind to GTEx
# system("echo 'Formatting combined data frame...'")
# tg <- tg[,-c(1:9), with=FALSE]
# d <- cbind(gtex,tg)
# # save as vcf
# system("echo 'Saving combined VCF file...'")
# vcf.total <- paste(outdir,"/rand.overlap.vcf",sep='')
# write.table(d, file=vcf.total, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
# # --------------------------------------------------------------------------------------------------------------

d <- fread(paste(outdir,"/rand.overlap.vcf",sep=''), sep='\t', header=TRUE, nrow=5)

## convert VCF to GDS
# system("echo 'Converting VCF to GDS...'")
# gds <- paste(outdir,"/rand.overlap.gds",sep='')
# snpgdsVCF2GDS(vcf.total, gds,  method="biallelic.only")
# # PCA the simple way 
# genofile <- snpgdsOpen(gds)
# system("echo 'Performing PCA...'")
# intersect_pca <- snpgdsPCA(genofile, num.thread=12)

# save(intersect_pca, file=paste(outdir,"/intersect-pca.RData",sep=''))
load(file=paste(outdir,"/intersect-pca.RData",sep=''))

# figure out population labels
#GTEx 
#RACE column:
#	1 = Asian
# 	2 = Black or African American
# 	3 = White
# 	4 = American Indian or Alaska Native
# 	99 = Unknown
# 	98 = Not Reported
system("echo 'Setting up population labels...'")
gpheno <- read.delim(file=phenofile, sep='\t', header=TRUE)
header <- as.data.frame(colnames(d)[-(1:9)], stringsAsFactors=FALSE)
colnames(header) <- "SUBJID"
header$race <- "1000 Genomes Ref"
header$label <- 0
gtexids <- grepl("GTEX",header$SUBJID)
header$label[gtexids] <- gpheno$RACE[match(header$SUBJID[gtexids],gpheno$SUBJID)]
header$race[gtexids] <- "OTH"
header$race[header$label == '1'] <- "ASN"
header$race[header$label == '2'] <- "AFR"
header$race[header$label == '3'] <- "EUR"
header$race[header$label == '4'] <- "IND"

#1000 Genomes
ref.codes <- read.delim(file=tg.pop.codes, sep='\t', stringsAsFactors=FALSE, header=TRUE) #c("Sample","Population")
ref.colors <- read.delim(file=tg.colors, sep='\t', stringsAsFactors=FALSE, header=TRUE) #c("hex","pop","superpop")
ref.codes$hex <- ref.colors$hex[match(ref.codes$Population, ref.colors$pop)]
#header$hex <- "#000000"
tgids <- !grepl("GTEX",header$SUBJID)
#header$hex[tgids] <- ref.codes$hex[match(header$SUBJID[tgids],ref.codes$Sample)]
header$pop <- "GTEx"
header$pop[tgids] <- ref.codes$Population[match(header$SUBJID[tgids],ref.codes$Sample)]
header$dataset <- "1000 Genomes"
header$dataset[grepl("GTEX",header$SUBJID)] <- "GTEx"
header$dataset <- as.factor(header$dataset)

# now plot
system("echo 'Plotting...'")

# get data frames that ggplot can work with
vcfdf <- as.data.frame(intersect_pca$eigenvect[,1:4])
colnames(vcfdf) <- c("PC1","PC2","PC3","PC4")
vcfdf$Race <- header$race
vcfdf$pop <- header$pop
vcfdf$dataset <- header$dataset
#vcfdf$hex <- as.factor(header$hex)
# reverse order of rows so that GTEx is plotted last
vcfdf <- vcfdf[rev(rownames(vcfdf)),]

# group.colors <- c(LWK="#CB9B31",
# 	GWD="#FBB916",
# 	MSL="#E0B719",
# 	ACB="#F3981C",
# 	ASW="#EB651E",
# 	YRI="#FCB839",
# 	ESN="#FFCD10",
# 	BEB="#811D80",
# 	STU="#A8529D",
# 	ITU="#B12F60",
# 	PJL="#DC1387",
# 	GIH="#6C3F94",
# 	CHB="#ABCA13",
# 	KHV="#4CAF38",
# 	CHS="#67B42C",
# 	JPT="#168D34",
# 	CDX="#369934",
# 	TSI="#272F6F",
# 	CEU="#234A9D",
# 	IBS="#678FCA",
# 	GBR="#66C8E2",
# 	FIN="#33BAC8",
# 	PEL="#E6131E",
# 	MXL="#DF0231",
# 	CLM="#CB3234",
# 	PUR="#BF3013",
# 	GTEx="#000000",
# 	black="black")

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

# to color by GTEx self-reported population
multi <- function(d,x,y){
	ggplot(data=d, aes(x=d[,x],y=d[,y],colour=Race)) + 
		geom_point(alpha=0.7, size=1) +
		scale_colour_manual(values = c("gray","purple","red","green","cyan","black"), name="Self-Reported Ancestry", 
			breaks=c("ASN","EUR","AFR","IND","OTH","1000 Genomes Ref")) +
		theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
  		xlab(paste("PC",x,sep='')) + ylab(paste("PC",y,sep=''))
}

jpeg("/users/nicolerg/gtex-admix/outputs/intersection-pca.jpg", width=2400, height=1600, res=100)
a <- multi(vcfdf,1,2)
b <- multi(vcfdf,1,3)
c <- multi(vcfdf,1,4)
d <- multi(vcfdf,2,3)
e <- multi(vcfdf,2,4)
f <- multi(vcfdf,3,4)
multiplot(a,b,c,d,e,f, cols=3)
dev.off()

a <- ggplot(data=vcfdf) +
	geom_point(data=subset(vcfdf, !pop=="GTEx"), aes(x=PC1, y=PC2, colour=pop, shape=dataset), alpha=1, size=1) + 
	theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
	scale_colour_manual(values=group.colors,name="1000 Genome Populations",breaks=breaks) +
	geom_point(data=subset(vcfdf, pop=="GTEx"), aes(x=PC1, y=PC2, shape=dataset), colour="black", size=1.5) +
	scale_shape_manual(values=c(20,2),name="Consortium") +
	xlab("PC1") + ylab("PC2")
b <- ggplot(data=vcfdf) +
	geom_point(data=subset(vcfdf, !pop=="GTEx"), aes(x=PC1, y=PC3, colour=pop, shape=dataset), alpha=1, size=1) + 
	theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
	scale_colour_manual(values=group.colors,name="1000 Genome Populations",breaks=breaks) +
	geom_point(data=subset(vcfdf, pop=="GTEx"), aes(x=PC1, y=PC3, shape=dataset), colour="black", size=1.5) +
	scale_shape_manual(values=c(20,2),name="Consortium") +
	xlab("PC1") + ylab("PC3")
c <- ggplot(data=vcfdf) +
	geom_point(data=subset(vcfdf, !pop=="GTEx"), aes(x=PC1, y=PC4, colour=pop, shape=dataset), alpha=1, size=1) + 
	theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
	scale_colour_manual(values=group.colors,name="1000 Genome Populations",breaks=breaks) +
	geom_point(data=subset(vcfdf, pop=="GTEx"), aes(x=PC1, y=PC4, shape=dataset), colour="black", size=1.5) +
	scale_shape_manual(values=c(20,2),name="Consortium") +
	xlab("PC1") + ylab("PC4")
d <- ggplot(data=vcfdf) +
	geom_point(data=subset(vcfdf, !pop=="GTEx"), aes(x=PC2, y=PC3, colour=pop, shape=dataset), alpha=1, size=1) + 
	theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
	scale_colour_manual(values=group.colors,name="1000 Genome Populations",breaks=breaks) +
	geom_point(data=subset(vcfdf, pop=="GTEx"), aes(x=PC2, y=PC3, shape=dataset), colour="black", size=1.5) +
	scale_shape_manual(values=c(20,2),name="Consortium") +
	xlab("PC2") + ylab("PC3")
e <- ggplot(data=vcfdf) +
	geom_point(data=subset(vcfdf, !pop=="GTEx"), aes(x=PC2, y=PC4, colour=pop, shape=dataset), alpha=1, size=1) + 
	theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
	scale_colour_manual(values=group.colors,name="1000 Genome Populations",breaks=breaks) +
	geom_point(data=subset(vcfdf, pop=="GTEx"), aes(x=PC2, y=PC4, shape=dataset), colour="black", size=1.5) +
	scale_shape_manual(values=c(20,2),name="Consortium") +
	xlab("PC2") + ylab("PC4")
f <- ggplot(data=vcfdf) +
	geom_point(data=subset(vcfdf, !pop=="GTEx"), aes(x=PC3, y=PC4, colour=pop, shape=dataset), alpha=1, size=1) + 
	theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
	scale_colour_manual(values=group.colors,name="1000 Genome Populations",breaks=breaks) +
	geom_point(data=subset(vcfdf, pop=="GTEx"), aes(x=PC3, y=PC4, shape=dataset), colour="black", size=1.5) +
	scale_shape_manual(values=c(20,2),name="Consortium") +
	xlab("PC3") + ylab("PC4")

jpeg("/users/nicolerg/gtex-admix/outputs/intersection-pca-refpop.jpg", width=1600, height=1200, res=200)
a
dev.off()

jpeg("/users/nicolerg/gtex-admix/outputs/intersection-pca-refpop-multiPC.jpg", width=4200, height=2000, res=200)
multiplot(a,b,c,d,e,f, cols=3)
dev.off()

system("echo 'Done.'")


# # FOR SEPARATE 1000 GENOME AND GTEX VCF FILES (NOT RECOMMENDED) ----------------------------------
# # first for 1000 Genomes

# gds <- paste(outdir,"/1000genomes.rand.all.gds",sep='')
# #snpgdsVCF2GDS(vcf.in.1kg, gds,  method="biallelic.only")
# #snpgdsSummary(gds)

# # PCA the simple way:
# #genofile <- openfn.gds(gds)
# #ccm_pca <- snpgdsPCA(genofile)
# #save(ccm_pca, file=paste(outdir,"/1kg-pca.RData",sep=''))
# load(paste(outdir,"/1kg-pca.RData",sep=''))


# # now for GTEx

# gds <- paste(outdir,"/gtex.rand.all.gds",sep='')
# #snpgdsVCF2GDS(vcf.in.gtex, gds,  method="biallelic.only")
# #snpgdsSummary(gds)
# #genofile <- openfn.gds(gds)
# #gtex_pca <- snpgdsPCA(genofile)
# #save(gtex_pca, file=paste(outdir,"/gtex-pca.RData",sep=''))
# load(paste(outdir,"/gtex-pca.RData",sep=''))

# # set things up for plotting

# # label GTEx by self-reported ancestry:
# ## Import phenotype tables
# # GTEx 
# # RACE column:
# #	1 = Asian
# # 	2 = Black or African American
# # 	3 = White
# # 	4 = American Indian or Alaska Native
# # 	99 = Unknown
# # 	98 = Not Reported
# gpheno <- read.delim(file=phenofile, sep='\t', header=TRUE)
# # get header from vcf
# con <- file(vcf.in.gtex,"r")
# header <- readLines(con,n=1)
# close(con)
# header <- data.frame(strsplit(header, '\t'), stringsAsFactors=FALSE)
# header <- header[-(1:9),]
# header <- data.frame(header, stringsAsFactors=FALSE)
# colnames(header) <- "SUBJID"
# header$label <- gpheno$RACE[match(header$SUBJID, gpheno$SUBJID)]
# header$race <- "OTH"
# header$race[header$label == '1'] <- "ASN"
# header$race[header$label == '2'] <- "AFR"
# header$race[header$label == '3'] <- "EUR"
# header$race[header$label == '4'] <- "IND"

# # make data frames that ggplot can work with
# gtexdf <- as.data.frame(gtex_pca$eigenvect[,1:4])
# colnames(gtexdf) <- c("PC1","PC2","PC3","PC4")
# gtexdf$Race <- header$race
# tgdf <- as.data.frame(ccm_pca$eigenvect[,1:4])
# colnames(tgdf) <- c("PC1","PC2","PC3","PC4")

# # first, plot just 1000 Genomes
# # with ggplot
# ggpc <- function(d,x,y){
# 	ggplot(data=d, aes(x=d[,x],y=d[,y])) +
# 		geom_point(alpha=0.7, size=1, colour="gray") +
# 		theme(plot.title=element_text(hjust = 0.5)) +
#   		labs(title=paste("PC", x, " vs. PC", y, sep='')) +
#   		xlab(paste("PC",x,sep='')) + ylab(paste("PC",y,sep=''))
# }

# # with base plot()
# plotpc <- function(a,b){
# 	plot(ccm_pca$eigenvect[,a],ccm_pca$eigenvect[,b],
# 		col="gray", 
# 		pch=20, 
# 		main=paste("PC", a, " vs. PC", b, sep=''), 
# 		xlab=paste("PC",a,sep=' '), 
# 		ylab=paste("PC",b,sep=' '))
# 	# points(gtex_pca$eigenvect[,a],gtex_pca$eigenvect[,b],
# 	# 	col="red", 
# 	# 	pch=20)
# }

# # first, plot just 1000 Genomes
# jpeg("/users/nicolerg/gtex-admix/outputs/1kg-pca.jpg", width=2600, height=1600, res=100)
# # par(mfrow=c(2,3))
# # plotpc(1,2)
# # plotpc(1,3)
# # plotpc(1,4)
# # plotpc(2,3)
# # plotpc(2,4)
# # plotpc(3,4)
# a <- ggpc(tgdf,1,2)
# b <- ggpc(tgdf,1,3)
# c <- ggpc(tgdf,1,4)
# d <- ggpc(tgdf,2,3)
# e <- ggpc(tgdf,2,4)
# f <- ggpc(tgdf,3,4)
# multiplot(a,b,c,d,e,f, cols=3)
# dev.off()

# # now, plot just GTEx
# # with ggplot
# ggtex <- function(d,x,y){
# 	ggplot(data=d, aes(x=d[,x],y=d[,y],colour=Race)) + 
# 		geom_point(alpha=0.7, size=1) +
# 		scale_colour_manual(values = c("purple","red","green","cyan","black")) +
# 		theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
#   		labs(title=paste("PC", x, " vs. PC", y, sep='')) +
#   		xlab(paste("PC",x,sep='')) + ylab(paste("PC",y,sep=''))
# }

# # with base plot()
# plotgtex <- function(a,b){
# 	plot(gtex_pca$eigenvect[,a],gtex_pca$eigenvect[,b],
# 		pch=20, 
# 		main=paste("PC", a, " vs. PC", b, sep=''), 
# 		xlab=paste("PC",a,sep=' '), 
# 		ylab=paste("PC",b,sep=' '),
# 		col=as.factor(header$race))
# 	# points(gtex_pca$eigenvect[,a],gtex_pca$eigenvect[,b],
# 	# 	col="red", 
# 	# 	pch=20)
# }

# jpeg("/users/nicolerg/gtex-admix/outputs/gtex-pca-2.jpg", width=2600, height=1600, res=100)
# # par(mfrow=c(2,3))
# # plotgtex(1,2)
# # plotgtex(1,3)
# # plotgtex(1,4)
# # plotgtex(2,3)
# # plotgtex(2,4)
# # plotgtex(3,4)
# # legend(0.4,0.6,unique(header$race), col=1:length(header$race), pch=20, cex=1.5)
# a <- ggtex(gtexdf,1,2)
# b <- ggtex(gtexdf,1,3)
# c <- ggtex(gtexdf,1,4)
# d <- ggtex(gtexdf,2,3)
# e <- ggtex(gtexdf,2,4)
# f <- ggtex(gtexdf,3,4)
# multiplot(a,b,c,d,e,f, cols=3)
# dev.off()

# # plot both with GTEx group labels:
# # with ggplot
# tgdf$Race <- "Ref"
# dat <- rbind(tgdf, gtexdf)

# plotboth <- function(comb,x,y){
# 	ggplot(data=comb, aes(x=comb[,x],y=comb[,y],colour=Race)) + 
# 		geom_point(alpha=0.7, size=1) +
# 		scale_colour_manual(values = c("purple","red","green","cyan","black","gray")) +
# 		theme(legend.position="right",plot.title=element_text(hjust = 0.5)) +
#   		labs(title=paste("PC", x, " vs. PC", y, sep='')) +
#   		xlab(paste("PC",x,sep='')) + ylab(paste("PC",y,sep=''))
# }

# # with base plot()
# plotpccol <- function(a,b){
# 	plot(ccm_pca$eigenvect[,a],ccm_pca$eigenvect[,b],
# 		col="gray", 
# 		pch=20, 
# 		main=paste("PC", a, " vs. PC", b, sep=''), 
# 		xlim=c(-0.03,0.1),
# 		ylim=c(-0.03,0.03),
# 		xlab=paste("PC",a,sep=' '), 
# 		ylab=paste("PC",b,sep=' '))
# 	points(gtex_pca$eigenvect[,a],gtex_pca$eigenvect[,b],
# 		col=as.factor(header$race), 
# 		pch=20)
# }

# jpeg("/users/nicolerg/gtex-admix/outputs/overlay-poplab-pca-2.jpg", width=2600, height=1600, res=100)
# # par(mfrow=c(2,3))
# # plotpccol(1,2)
# # plotpccol(1,3)
# # plotpccol(1,4)
# # plotpccol(2,3)
# # plotpccol(2,4)
# # plotpccol(3,4)
# # legend(0.08,0.0,unique(header$race), col=1:length(header$race), pch=20, cex=2)
# a <- plotboth(dat,1,2)
# b <- plotboth(dat,1,3)
# c <- plotboth(dat,1,4)
# d <- plotboth(dat,2,3)
# e <- plotboth(dat,2,4)
# f <- plotboth(dat,3,4)
# multiplot(a,b,c,d,e,f, cols=3)
# dev.off()

# # # a more thorough way (source: http://corearray.sourceforge.net/tutorials/SNPRelate/#data-formats-used-in-snprelate)
# # # make pop_groups file (population abbreviation for every sample)
# # # import the two pop_codes tables
# # sample_pop <- read.delim(file=tg.pop.codes, sep='\t', header=TRUE)
# # pop_match <- read.delim(file=tg.pop.match, sep='\t', header=TRUE)
# # # import the vcf header to get order of samples
# # con <- file(vcf.in.1kg,"r")
# # header <- readLines(con,n=1)
# # close(con)
# # header <- data.frame(strsplit(header, '\t'), stringsAsFactors=FALSE)
# # header <- header[-(1:9),]
# # header <- data.frame(header, stringsAsFactors=FALSE)
# # colnames(header) <- "sample_id"

# # pop_code <- read.gdsn(index.gdsn(genofile, path=pop_groups))

# # table(pop_code)
# # head(pop_code)
# # set.seed(1000)
# # snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
# # snpset.id <- unlist(snpset)

# # #PCA
# # pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=12)
# # pc.percent <- pca$varprop*100
# # head(round(pc.percent, 2))
# # # Get sample id
# # sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# # # assume the order of sample IDs is as the same as population codes
# # head(cbind(sample.id, pop_code))

# # tab <- data.frame(sample.id = pca$sample.id,
# #     pop = factor(pop_code)[match(pca$sample.id, sample.id)],
# #     EV1 = pca$eigenvect[,1],    # the first eigenvector
# #     EV2 = pca$eigenvect[,2],    # the second eigenvector
# #     stringsAsFactors = FALSE)
# # head(tab)

# # plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
# # legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

# # # Plot the principal component pairs for the first four PCs
# # lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
# # pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)