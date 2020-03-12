library(ggplot2)
library(reshape)

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

#load("/mnt/lab_data/montgomery/nicolerg/pca-subset/intersect-pca.RData") #intersect_pca
pcs <- read.delim("/mnt/lab_data/montgomery/nicolerg/pca-subset/gtex_PCs_from_intersect.csv",header=TRUE,stringsAsFactors=FALSE,sep=',')
globalai <- read.delim("/mnt/lab_data/montgomery/nicolerg/admix/gtex_ai_global.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)

pcs <- pcs[match(globalai$ID,rownames(pcs)),]

rownames(globalai) <- globalai$ID
globalai$ID <- NULL

df <- cbind(globalai,pcs)

# plotPCs <- function(pop,low,high,pc){
# 	ggplot(df,aes_string(x=pop,y=pc,colour=pop)) + 
# 		geom_point(size=0.8) + 
# 		scale_colour_gradient(low = low, high = high) +
# 		labs(x=paste("Fraction",pop,"Global Ancestry",sep=' '))
# }
head(df)
ddf <- melt(df, id=c("ASN","EUR","AFR","UNK"))
head(ddf)
colnames(ddf) <- c("ASN","EUR","AFR","UNK","pc_var","pc_val")

plotPCs <- function(pop,low,high){
	ggplot(ddf,aes_string(x=pop,y="pc_val",colour=pop)) + 
		geom_point(size=0.8) + 
		scale_colour_gradient(low = low, high = high) +
		labs(x=paste("Fraction",pop,"Global Ancestry",sep=' ')) +
		facet_wrap("pc_var", scales="free", ncol=4)
}

a <- plotPCs("ASN","#bdff5b","#456614")
b <- plotPCs("AFR","#ffcc7f","#7f4c00")
c <- plotPCs("EUR","#b1f5ff","#285e66")

jpeg("/users/nicolerg/gtex-admix/pcs_vs_globalai_facet.jpg", width=2600, height=2000, res=200)
multiplot(a,b,c, cols=1)
dev.off()

#### not facet-wrapped
# #ASN
# a <- plotPCs("ASN","#bdff5b","#456614","PC1")
# b <- plotPCs("ASN","#bdff5b","#456614","PC2")
# c <- plotPCs("ASN","#bdff5b","#456614","PC3")
# #AFR
# d <- plotPCs("AFR","#ffcc7f","#7f4c00","PC1")
# e <- plotPCs("AFR","#ffcc7f","#7f4c00","PC2")
# f <- plotPCs("AFR","#ffcc7f","#7f4c00","PC3")
# #EUR
# g <- plotPCs("EUR","#b1f5ff","#285e66","PC1")
# h <- plotPCs("EUR","#b1f5ff","#285e66","PC2")
# i <- plotPCs("EUR","#b1f5ff","#285e66","PC3")

# jpeg("/users/nicolerg/gtex-admix/pcs_vs_globalai.jpg", width=3000, height=2000, res=200)
# multiplot(a,b,c,d,e,f,g,h,i, cols=3)
# dev.off()
####

# plot by population
pdata <- melt(df, id=c("ASN","EUR","AFR","UNK"))
jpeg("/users/nicolerg/gtex-admix/pcs_vs_asn.jpg", width=1600, height=1200, res=200)
ggplot(pdata,aes(x=ASN,y=value,colour=ASN)) + 
	geom_point(size=0.8) + 
	scale_colour_gradient(low = "#bdff5b", high = "#456614") +
	labs(x=paste("Fraction Asian Global Ancestry",sep=' '),y="") +
	facet_wrap(~variable,scales="free")
dev.off()

jpeg("/users/nicolerg/gtex-admix/pcs_vs_eur.jpg", width=1600, height=1200, res=200)
ggplot(pdata,aes(x=EUR,y=value,colour=EUR)) + 
	geom_point(size=0.8) + 
	scale_colour_gradient(low = "#b1f5ff", high = "#285e66") +
	labs(x=paste("Fraction European Global Ancestry",sep=' '),y="") +
	facet_wrap(~variable,scales="free")
dev.off()

jpeg("/users/nicolerg/gtex-admix/pcs_vs_afr.jpg", width=1600, height=1200, res=200)
ggplot(pdata,aes(x=AFR,y=value,colour=AFR)) + 
	geom_point(size=0.8) + 
	scale_colour_gradient(low = "#ffcc7f", high = "#7f4c00") +
	labs(x=paste("Fraction African Global Ancestry",sep=' '),y="") +
	facet_wrap(~variable,scales="free")
dev.off()

#### not facet-wrapped
## plot all groups on one plot
# mdata <- melt(df, id=c("PC1","PC2","PC3","PC4"))

# a <- ggplot(mdata,aes(x=value,y=PC1)) + 
# 	geom_point(aes(fill=variable,size=value),colour="black",alpha=0.8, pch=21) +
# 	scale_fill_manual(values=c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="#000000"),name="Local Ancestry") +
# 	scale_size_continuous(range = c(0.1, 3)) +
# 	labs(x="Fraction of Global Ancestry")
# b <- ggplot(mdata,aes(x=value,y=PC2)) + 
# 	geom_point(aes(fill=variable,size=value),colour="black",alpha=0.8, pch=21) +
# 	scale_fill_manual(values=c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="#000000"),name="Local Ancestry") +
# 	scale_size_continuous(range = c(0.1, 3)) +
# 	labs(x="Fraction of Global Ancestry")
# c <- ggplot(mdata,aes(x=value,y=PC3)) + 
# 	geom_point(aes(fill=variable,size=value),colour="black",alpha=0.8, pch=21) +
# 	scale_fill_manual(values=c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="#000000"),name="Local Ancestry") +
# 	scale_size_continuous(range = c(0.1, 3)) +
# 	labs(x="Fraction of Global Ancestry")
# d <- ggplot(mdata,aes(x=value,y=PC4)) + 
# 	geom_point(aes(fill=variable,size=value),colour="black",alpha=0.8, pch=21) +
# 	scale_fill_manual(values=c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="#000000"),name="Local Ancestry") +
# 	scale_size_continuous(range = c(0.1, 3)) +
# 	labs(x="Fraction of Global Ancestry")

# jpeg("/users/nicolerg/gtex-admix/pcs_vs_allglobalai.jpg", width=2400, height=1800, res=200)
# multiplot(a,b,c,d,cols=2)
# dev.off()
####

# plot all groups on one plot
mdata <- melt(df, id=c("ASN","EUR","AFR","UNK"))
head(mdata)
colnames(mdata) <- c("ASN","EUR","AFR","UNK","PCs","PC_val")
mdata <- melt(mdata, id=c("PCs","PC_val"))
head(mdata)
mdata <- mdata[!grepl("UNK",mdata$variable),]

jpeg("/users/nicolerg/gtex-admix/pcs_vs_allglobalai_facet.jpg", width=2200, height=1800, res=200)
ggplot(mdata,aes(x=value,y=PC_val)) + 
	geom_point(aes(fill=variable,size=value),colour="black",alpha=0.8, pch=21) +
	scale_fill_manual(values=c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="#000000"),name="Local Ancestry") +
	scale_size_continuous(range = c(0.1, 3)) +
	labs(x="Fraction of Global Ancestry") +
	facet_wrap(~PCs,scales="free")
dev.off()
