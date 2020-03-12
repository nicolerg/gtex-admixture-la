library(ggplot2)

# read in global ancestry proportions
indir <- "/mnt/lab_data/montgomery/nicolerg/rfmix-admixed/combined/"
global_ai <- read.delim(paste0(indir,"gtex_ai_global_all.txt"), sep='\t', header=TRUE, stringsAsFactors=FALSE)

# filter out subject with >0.9 homogeneity
global_ai <- global_ai[ global_ai$AFR < 0.9 & global_ai$ASN < 0.9 & global_ai$EUR < 0.9 ,]
admixed <- data.frame(ID=global_ai$ID)
write.table(admixed,paste0(indir,"gtex-admixed0.9.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
n <- nrow(admixed)
# ID	ASN	EUR	AFR
# GTEX-1117F	0.0064	0.1938	0.7997

master <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(master) <- c("chr","start","stop","pop","ind","eur")

inds <- global_ai$ID

for (ind in inds){
	bedA <- paste0(indir,ind,"_A.bed")
	bedB <- paste0(indir,ind,"_B.bed")

	labs <- c("Chr","Start","Stop","Pop","Ind")
	df1 <- read.delim(bedA, sep='\t', header=FALSE, stringsAsFactors=FALSE)
	df1 <- df1[,1:4]
	df1$Ind <- paste(ind,"_A",sep='')
	colnames(df1) <- labs
	df2 <- read.delim(bedB, sep='\t', header=FALSE, stringsAsFactors=FALSE)
	df2 <- df2[,1:4]
	df2$Ind <- paste(ind,"_B",sep='')
	colnames(df2) <- labs

	df <- data.frame(rbind(df1, df2), stringsAsFactors=FALSE)
	df$Eur <- global_ai$EUR[match(ind, global_ai$ID)]

	# now add to master df
	master <- data.frame(rbind(master, df), stringsAsFactors=FALSE)
}

# standardize start positions and make a table of max stop values for each chromosome
dict <- data.frame(c(1:22),rep(NA,22),rep(NA,22),rep(NA,22))
colnames(dict) <- c("Chr","Max","Start","Stop")
for (i in c(1:22)){
	rows <- which(master$Chr == i)
	min <- min(master$Start[rows])
	master$Start[rows] <- master$Start[rows] - min # make each chromosome start at 0
	master$Stop[rows] <- master$Stop[rows] - min # shift all segments over by that amount
	max <- max(master$Stop[rows])
	dict$Max[i] <- max
}

# figure out start and stop positions for each chromosome
for (chr in c(1:22)){
	if (chr == 1){
		dict$Start[chr] <- 0
		dict$Stop[chr] <- dict$Max[chr]
	} else {
		dict$Start[chr] <- dict$Stop[chr-1] # make chromosomes start and stop on the same point
		dict$Stop[chr] <- dict$Start[chr] + dict$Max[chr]
	}
}

# add start value to every "start" and "stop" in master
for (chr in c(1:22)){
	rows <- which(master$Chr == chr)
	master$Start[rows] <- master$Start[rows] + dict$Start[chr]
	master$Stop[rows] <- master$Stop[rows] + dict$Start[chr]
}

# sort master by descending eur values
master <- master[order(master$Eur, decreasing=TRUE),] 
# assign a y value to each haplotype (master$ind)
master$Y <- NA
indiv <- unique(master$Ind) #IND_A, #IND_B
yval <- data.frame(indiv,c(1:length(indiv)),stringsAsFactors=FALSE)
colnames(yval) <- c("ind","y")
matchID <- function(x){
	yval$y[grep(x,yval$ind)]
}
master$Y <- unlist(lapply(master$Ind, matchID))
master$Y <- as.numeric(master$Y)

get_id <- function(x){
	return(unlist(unname(strsplit(x,'_')))[1])
}
ids <- unique(sapply(master$Ind, get_id))
breaks <- seq(1.5,233.5,by=2)
yint <- seq(2.5,232.5,by=2)

save(master, dict, file="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/figure_data/block-admix-haplo.RData")

pdf("/users/nicolerg/gtex-admix/plots/final_figures/block-admix-haplo.pdf", width=9, height=8)
ggplot(master) + 
	geom_segment(aes(x=Start, xend=Stop, y=Y, yend=Y, colour=Pop), size=1, lineend="butt") +
	geom_hline(yintercept=yint, size=0.5) +
	geom_point(x=NA,y=NA,shape=22,aes(fill=Pop),colour='black',size=5) +
	theme_bw() +
	theme(panel.border = element_rect(size=1, colour='black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          #axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 4, colour='black'),
          axis.text.x = element_text(size = 8, colour='black'),
          plot.title = element_text(hjust = 0.5, size = 14, colour='black'),
          legend.text = element_text(size = 12, colour='black'),
          legend.title = element_blank(),
          legend.position = 'bottom',
          axis.title.x = element_text(size = 14, colour='black'),
          axis.title.y = element_blank()) +
	scale_colour_manual(name="Local Ancestry", 
		values = c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
		guide = 'none') +
	scale_fill_manual(name="Local Ancestry", 
		values = c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
		breaks=c("ASN","EUR","AFR","UNK"),
		labels=c("Asian","European","African","Unknown")) +
	geom_vline(xintercept=dict$Stop[1:21], size=0.5) + 
	scale_x_discrete(limits=(dict$Start+dict$Stop)/2,
        labels=c(1:22),name="Chromosome") +
	scale_y_continuous(breaks=breaks, labels=ids, expand = c(0,0), limits=c(0.5,234.5)) +
  	#scale_y_continuous(expand = c(0, 0)) +
  	labs(title=paste0("GTEx v8 Admixed Individuals (N=",n,")"),y=NULL) 
dev.off()

# # for paper figure

# df <- data.frame(c(rep("White",714),rep("African_American",103),rep("Asian",11)))
# colnames(df) <- "Race"
# df$y <- c(1:nrow(df))
# df$xstart <- 0
# df$xend <- 1

# jpeg("/users/nicolerg/gtex-admix/plots/self-reported-block.jpg", width=200, height=700, res=200)
# ggplot(df) + geom_segment(aes(x=xstart, xend=xend, y=y, yend=y, colour=Race), size=0.8, lineend="butt") +
# 	theme_bw() + 
# 	theme(panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), 
#           axis.ticks = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_text(size = 14)) +
# 	scale_colour_manual(values = c(Asian="#ADFF33",White="#0000FF",African_American="#FF9900")) + 
# 	labs(x = NULL, y = "828 Individuals")
# dev.off()
