# written for desktop. make heatmap of correlation between PCs and PEER factors

setwd("~/Documents/Surya/gtex-admix")
dat <- read.delim("Adipose_Subcutaneous.v8.covariates.txt",sep='\t', header=TRUE, stringsAsFactors = FALSE)
rownames(dat) <- dat$ID
dat$ID <- NULL

dat <- data.frame(t(dat),stringsAsFactors = FALSE)

lm.fit <- lm(PC1 ~. -PC2-PC3-PC4-PC5, data=dat)
summary(lm.fit)

attach(dat)
correlationPC1 <- rep(NA, 60)
for(i in 1:60){
  correlationPC1[i] <- cor(PC1, dat[,5+i])
}
correlationPC2 <- rep(NA, 60)
for(i in 1:60){
  correlationPC2[i] <- cor(PC2, dat[,5+i])
}
correlationPC3 <- rep(NA, 60)
for(i in 1:60){
  correlationPC3[i] <- cor(PC3, dat[,5+i])
}
correlationPC4 <- rep(NA, 60)
for(i in 1:60){
  correlationPC4[i] <- cor(PC4, dat[,5+i])
}
correlationPC5 <- rep(NA, 60)
for(i in 1:60){
  correlationPC5[i] <- cor(PC5, dat[,5+i])
}
corr <- data.frame(correlationPC1,correlationPC2,correlationPC3,correlationPC4,correlationPC5)
colnames(corr) <- c("PC1","PC2","PC3","PC4","PC5")

#heatmap
library(gplots)
corrm <- as.matrix(corr)

jpeg('pc-peer-corr-adipose.jpg', width=600, height=800, res=100)
heatmap.2(corrm, trace='none', Rowv=FALSE, Colv=FALSE, dendrogram ='none', scale='none', col=rev(redblue(32)), keysize=1.1, breaks=seq(-1, 1, length.out=33), ylab='PEER Factor', key.title=NA)
dev.off()

# plot(corr$correlationPC1, type='l', col="green", xlab="PEER Factor Index", ylab="Correlation with Genotype PC", ylim=c(-0.4,.4))
# lines(corr$correlationPC2, col="blue")
# lines(corr$correlationPC3, col="red")
# lines(corr$correlationPC4, col="purple")
# lines(corr$correlationPC5, col="orange")
# legend("topright", c('PC1','PC2', 'PC3', 'PC4', 'PC5'), col=c('green','blue','red','purple','orange'), lty=1, cex=0.6)