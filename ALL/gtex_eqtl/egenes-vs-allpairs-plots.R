library(ggplot2)

##### first, make plots for just eGene comparisons

mypairs <- read.delim("/users/nicolerg/gtex-admix/egenes_overlap_v8.tsv", sep='\t', header=TRUE, stringsAsFactors=FALSE)

# scatterplot of pvals
sub <- mypairs[complete.cases(mypairs),]
sub$mylog <- -log(sub$my_pval_nominal)
sub$gtexlog <- -log(sub$gtex_pval)

fit <- lm(gtexlog ~ mylog, data=sub)
r <- summary(fit)$r.squared

print(r)

jpeg("/users/nicolerg/gtex-admix/pval-correlation.jpg", width=1600, height=1400, res=200)
ggplot(sub, aes(x=mylog,y=gtexlog)) + 
	geom_point() + 
	annotate("text", x = 100, y = 550, label = paste("R^2:",round(r,digits=3),sep=' ')) +
	labs(x="-log(my_nominal_pval)",y="-log(GTEx_nominal_pval)")
dev.off()

# histogram of pvals not included in GTex
sub <- mypairs[!complete.cases(mypairs),]
sub$mylog <- -log(sub$my_pval_nominal)
jpeg("/users/nicolerg/gtex-admix/my-pval-hist.jpg", width=1600, height=1400, res=200)
hist(as.numeric(sub$mylog),breaks=2000,xlab="-log(my_pval_nominal)",xlim=c(0,60),main="Unshared pvals (my version)")
dev.off()

# histogram of GTEx pvals
gtexpairs_file="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz"
gtexpairs <- read.delim(gtexpairs_file, sep='\t', header=TRUE, stringsAsFactors=FALSE)
jpeg("/users/nicolerg/gtex-admix/gtex-pval-hist.jpg", width=1600, height=1400, res=200)
hist(as.numeric(-log(gtexpairs$pval_nominal)),breaks=2000,xlab="-log(GTEx_pval_nominal)",main="All v8 GTEx pvals",xlim=c(0,60))
dev.off()

##### now make plots for my eGenes vs GTEx all pairs file

pval_file <- "/users/nicolerg/gtex-admix/allpairs_overlap_v8.tsv"
pvals <- read.delim(pval_file, sep='\t', header=TRUE, stringsAsFactors=FALSE)

# scatterplot of pvals
pvals$mylog <- -log(pvals$my_pval_nominal)
pvals$gtexlog <- -log(pvals$gtex_pval_nominal)

fit <- lm(gtexlog ~ mylog, data=pvals)
r <- summary(fit)$r.squared

print(r)

jpeg("/users/nicolerg/gtex-admix/pval-correlation-unshared.jpg", width=1600, height=1400, res=200)
ggplot(pvals, aes(x=mylog,y=gtexlog)) + 
	geom_point() + 
	annotate("text", x = 100, y = 550, label = paste("R^2:",round(r,digits=3),sep=' ')) +
	labs(x="-log(my_nominal_pval)",y="-log(GTEx_nominal_pval)")
dev.off()

# instead, make histogram of all gtex_pval from all_pairs 
jpeg("/users/nicolerg/gtex-admix/gtex-overlap-pval-hist.jpg", width=1600, height=1400, res=200)
hist(as.numeric(pvals$gtexlog),breaks=2000,xlab="-log(gtex_pval_nominal)",xlim=c(0,60),main="GTEx All Pair p values for my matching eGenes")
dev.off()