library(data.table)
library(ggplot2)
library(OneR)

# 29 April 2019 
# see how much the number of SNPs tested for a gene predicts whether or not there is a different lead SNP

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-06.RData')

egenes_master <- egenes_master[pval_nominal < 1e-6]
diff_lead <- master_diff_lead

egenes_master[,difflead := ifelse(gene_id %in% diff_lead[,gene_id], 1, 0)]
sub <- egenes_master[,.(gene_id, count, difflead, tissue)]
sub <- unique(sub)
glm.fit <- glm(difflead ~ I(count^2) + count + tissue, data=sub, family='binomial')
summary(glm.fit)

logit.plot = function(data,predictor,outcome,numIntervals=10,spline=TRUE,spar=0.7) {
  #Check inputs
  if (!is.character(predictor)) {stop("The predictor must be a string.")}
  if (!is.character(outcome)) {stop("The outcome must be a string.")}
  
  #Initialize vectors
  logits = rep(0,numIntervals)
  means = rep(0,numIntervals)
  
  #Define outcome and predictor columns (renames them to make the code below simpler)
  outcomeIndex = which(colnames(data)==outcome)
  predictorIndex = which(colnames(data)==predictor)
  names(data)[outcomeIndex] = "outcome"
  names(data)[predictorIndex] = "predictor"
  
  #Check for missing data
  if (any(is.na(data$predictor))) {stop("The predictor contains missing data.")}
  if (any(is.na(data$outcome))) {stop("The outcome contains missing data.")}
  
  #Sort data in ascending order of predictor
  data=data[order(data$predictor),]
  
  #Define the intervals
  intervalLength = floor(dim(data)[1] / numIntervals)
  #floor() enforces rounding down, in case numIntervals doesn't divide evenly into dim(data)[1]
  
  if ((dim(data)[1] / numIntervals) - floor(dim(data)[1] / numIntervals) != 0) {
    warning("The number of intervals does not divide evenly into the size of the dataset. Some data points will be omitted.")
  }
  
  #Define the starting index of each section
  intervalStarts = c()
  for (k in 1:numIntervals) {
    intervalStarts[k] = 1 + (k-1)*intervalLength
  }
  
  #Loop over each section
  for (j in intervalStarts) {
    positive=0
    negative=0
    sum=0
    
    #Loop over each data point in the section
    for (i in c(j:(j+intervalLength-1))) {
      if (data$outcome[i] == 1) {             #outcome variable is 1
        positive = positive + 1
      } else if (data$outcome[i] == 0) {      #outcome variable is 0
        negative = negative + 1
      } else {
        stop("The outcome column must be binary.")
      }
      
      sum = sum + data$predictor[i]   #adding the predictor values so we can calculate the mean
    }
    
    x = positive/(positive+negative)        #calculates the probability
    logits[1+round(j/intervalLength)] = log(x/(1-x))        #puts the logit in the index of j
    means[1+round(j/intervalLength)] = sum/intervalLength   #puts the mean in the index of j 
    
  }
  
  #Check for infinite logits
  if (any(is.infinite(logits))) {
    warning("Infinite logits generated. Not all points will be plotted.")
  }
  
  #Plot each logit at the corresponding predictor value
  plot(means,logits,xlab=paste(predictor))
  
  #Plot a smoothing spline
  if (spline==TRUE) {
    spline = smooth.spline(means[!is.infinite(logits)], logits[!is.infinite(logits)], spar=spar)
    lines(spline)  
  }
}

pdf('~/gtex-admix/plots/final_figures/counts-vs-difflead.pdf',width=4,height=4)
logit.plot(as.data.frame(sub),'count','difflead',numIntervals=25,spline=TRUE,spar=0.7) 
dev.off()

# cyp2e1
egenes_master[grepl('ENSG00000130649',gene_id)]
# hla-drb6
egenes_master[grepl('ENSG00000229391',gene_id)]

# HLA-DRB6
global <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/input/skin_global_HLA.txt'
local <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/input/skin_local_HLA.txt'
self_reported_asthma <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/gtex_gwas/coloc_input/UKB_20002_1111_self_reported_asthma_HLA.txt'
diagnosed_asthma <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/gtex_gwas/coloc_input/UKB_6152_8_diagnosed_by_doctor_Asthma_HLA.txt'
tagc_asthma <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/gtex_gwas/coloc_input/TAGC_Asthma_EUR_HLA.txt'

global <- fread(global, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','gene'))
local <- fread(local, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','gene'))
self_reported_asthma <- fread(self_reported_asthma, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','effect_af'))
diagnosed_asthma <- fread(diagnosed_asthma, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','effect_af'))
tagc_asthma <- fread(tagc_asthma, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','effect_af'))

lead_snps <- data.frame(snp_pos = c(32660781,32660782,32660783,32604534), pvalue=c(1.723884e-07,1.723884e-07,1.723884e-07,9.305206e-08), method=c('Local','Local','Local','Global'))

pdf('~/gtex-admix/plots/final_figures/HLA-DRB6-coloc-locus.pdf',width=8,height=3)
ggplot() +
	geom_point(data=self_reported_asthma, aes(x=snp_pos, y=-log10(pvalue), colour='self_reported_asthma'), size=0.5, alpha=0.5) +
	geom_point(data=diagnosed_asthma, aes(x=snp_pos, y=-log10(pvalue), colour='diagnosed_asthma'), size=0.5, alpha=0.5) +
	geom_point(data=tagc_asthma, aes(x=snp_pos, y=-log10(pvalue), colour='tagc_asthma'), size=0.5, alpha=0.5) +
	geom_point(data=global, aes(x=snp_pos, y=-log10(pvalue), colour='Global'), size=0.5, alpha=0.5) +
	geom_point(data=local, aes(x=snp_pos, y=-log10(pvalue), colour='Local'), size=0.5, alpha=0.5) +
	geom_point(data=lead_snps, aes(x=snp_pos, y=-log10(pvalue), fill=factor(method)), colour='black', alpha=1, size=1.5, show.legend=FALSE, shape=21) +
	labs(x='Position on chr6',y='-log10(p-value)',title='HLA-DRB6 colocalization locus') +
	theme_bw() + 
	scale_colour_manual(values=c(Global='#8DD3C7',
								Local='#FDB462',
								self_reported_asthma='green3',
								diagnosed_asthma='purple',
								tagc_asthma='red')) +
	theme(legend.title=element_blank())
dev.off()


# CYP2E1
global <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/input/skin_global_CYP2E1.txt'
local <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/input/skin_local_CYP2E1.txt'
bmi <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/gtex_gwas/coloc_input/UKB_21001_Body_mass_index_BMI_CYP2E1.txt'
bfp <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/gtex_gwas/coloc_input/UKB_23099_Body_fat_percentage_CYP2E1.txt'

global <- fread(global, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','gene'))
local <- fread(local, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','gene'))
bmi <- fread(bmi, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','gene'))
bfp <- fread(bfp, sep='\t', header=FALSE, col.names=c('chr','snp_pos','alt','ref','beta','se','pvalue','gene'))

lead_snps <- data.frame(snp_pos = c(132171412,131077337), pvalue=c(2.00943e-07, 5.16194e-07), method=c('Local','Global'))

pdf('~/gtex-admix/plots/final_figures/CYP2E1-coloc-locus.pdf',width=8,height=3)
ggplot() +
	geom_point(data=bmi, aes(x=snp_pos, y=-log10(pvalue), colour='bmi'), size=0.5, alpha=0.5) +
	geom_point(data=bfp, aes(x=snp_pos, y=-log10(pvalue), colour='bfp'), size=0.5, alpha=0.5) +
	geom_point(data=global, aes(x=snp_pos, y=-log10(pvalue), colour='Global'), size=0.5, alpha=0.5) +
	geom_point(data=local, aes(x=snp_pos, y=-log10(pvalue), colour='Local'), size=0.5, alpha=0.5) +
	geom_point(data=lead_snps, aes(x=snp_pos, y=-log10(pvalue), fill=factor(method)), colour='black', alpha=1, size=1.5, show.legend=FALSE, shape=21) +
	labs(x='Position on chr10',y='-log10(p-value)',title='CYP2E1 colocalization locus') +
	theme_bw() + 
	scale_colour_manual(values=c(Global='#8DD3C7',
								Local='#FDB462',
								bmi='green3',
								bfp='purple')) +
	theme(legend.title=element_blank())
dev.off()

