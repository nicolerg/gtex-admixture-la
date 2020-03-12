# Author: Nicole Gay
# Winter 2018
# Plot outputs from slow_eqtl.R --> match_gtex_manual.py (nominal P-value scatterplots)

library(ggplot2)
library(data.table)

# args <- commandArgs(trailingOnly=TRUE)
# gtex_local <- args[1]
# scatter <- args[2]
# outfile <- args[3]

gtex_local="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/lung-chr1-gtex-local.tsv.gz"
scatter <- "/users/nicolerg/gtex-admix/lung/local-v-gtex-lung-chr1.jpg"
scatter_sub <- "/users/nicolerg/gtex-admix/lung/local-v-gtex-lung-chr1-zoom.jpg"
outfile <- "/users/nicolerg/gtex-admix/lung/lung-chr1-sigdiff-pairs.tsv"

# global vs local

cat("Loading data\n")
pvals <- fread(input = sprintf("zcat %s", gtex_local), sep='\t', header=TRUE, stringsAsFactors=FALSE)
pvals <- data.frame(pvals, stringsAsFactors=FALSE)

cat("Transforming p-values\n")
pvals <- pvals[complete.cases(pvals),]
pvals$local <- -log10(pvals$local_pval_nominal)
pvals$global <- -log10(pvals$gtex_pval_nominal)

# cat("Performing regression\n")
# lmfit <- lm(local ~ global, data=pvals)
# out <- summary(lmfit)$coefficients

# cat("Extracting coefficients\n")
# coef <- as.numeric(out[2][1])
# intc <- as.numeric(out[1][1])

bon_thresh <- -log10(0.05/nrow(pvals))
print(bon_thresh)

# cat("Defining outliers\n")
# pvals$upper <- pvals$global * coef + (intc + 3.8)
# pvals$lower <- pvals$global * coef - (intc + 3.8)
# pvals$outlier <- 1
# head(pvals)
# bin <- (pvals$local < pvals$lower) | (pvals$local > pvals$upper)
# pvals$outlier[bin] <- 2
# pvals$outlier[pvals$local < bon_thresh] <- 1
# pvals$outlier[pvals$global < bon_thresh] <- 1

# table(pvals$outlier)

cat("Making scatterplot\n")
jpeg(scatter, width=1000, height=800, res=150)
ggplot(pvals, aes(x=global,y=local)) + 
	geom_point(size=1.0) + 
	labs(x="Global Nominal P-Val (-log10)",y="Local Nominal P-Val (-log10)") +
	theme(legend.position="none") +
	geom_vline(xintercept=bon_thresh,linetype="twodash",size=1) +
	geom_hline(yintercept=bon_thresh,linetype="twodash",size=1) 
dev.off()

jpeg(scatter_sub, width=1000, height=800, res=150)
ggplot(pvals, aes(x=global,y=local)) + 
	geom_point(size=1.0) + 
	labs(x="Global Nominal P-Val (-log10)",y="Local Nominal P-Val (-log10)") +
	ylim(c(0,50)) +
	xlim(c(0,50)) +
	theme(legend.position="none") +
	geom_vline(xintercept=bon_thresh,linetype="twodash",size=1) +
	geom_hline(yintercept=bon_thresh,linetype="twodash",size=1) 
dev.off()

# jpeg(scatter, width=1000, height=800, res=150)
# ggplot(pvals, aes(x=global,y=local, colour=factor(outlier))) + 
# 	geom_point(size=1.0) + 
# 	labs(x="Global Nominal P-Val (-log10)",y="Local Nominal P-Val (-log10)") +
# 	scale_colour_manual(values=c("black","red")) +
# 	theme(legend.position="none") +
# 	geom_vline(xintercept=bon_thresh,linetype="twodash",size=1) +
# 	geom_hline(yintercept=bon_thresh,linetype="twodash",size=1) 
# dev.off()

# jpeg(scatter_sub, width=1000, height=800, res=150)
# ggplot(pvals, aes(x=global,y=local, colour=factor(outlier))) + 
# 	geom_point(size=1.0) + 
# 	labs(x="Global Nominal P-Val (-log10)",y="Local Nominal P-Val (-log10)") +
# 	scale_colour_manual(values=c("black","red")) +
# 	ylim(c(0,50)) +
# 	xlim(c(0,50)) +
# 	theme(legend.position="none") +
# 	geom_vline(xintercept=bon_thresh,linetype="twodash",size=1) +
# 	geom_hline(yintercept=bon_thresh,linetype="twodash",size=1) 
# dev.off()

# pvals <- pvals[pvals$outlier==2,1:2]
# write.table(pvals,file=outfile,sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)

## global vs gtex (baseline)

# pvals2 <- read.delim(global_gtex, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# # scatterplot of pvals
# pvals2$global <- -log(pvals2$global_pval_nominal)
# pvals2$gtex <- -log(pvals2$gtex_pval_nominal)

# fit <- lm(global ~ gtex, data=pvals2)
# r <- summary(fit)$r.squared

# print(r)

# jpeg("/users/nicolerg/gtex-admix/global-v-gtex-wholeblood-chr21.jpg", width=1600, height=1400, res=200)
# ggplot(pvals2, aes(x=global,y=gtex)) + 
# 	geom_point() + 
# 	annotate("text", x = 100, y = 250, label = paste("R-Squared:",round(r,digits=3),sep=' ')) +
# 	labs(x="Global Nominal P-Val (-log10)",y="GTEx Nominal P-Val (-log10)")
# dev.off()

# jpeg("/users/nicolerg/gtex-admix/global-v-gtex-hg19-500-noPC-zoom.jpg", width=1600, height=1400, res=200)
# ggplot(pvals2, aes(x=gtex,y=global)) + 
# 	geom_point() + 
# 	labs(x="-log(gtex_nominal_pval)",y="-log(global_nominal_pval)") +
# 	ylim(c(0,30)) +
# 	xlim(c(0,30))
# dev.off()
