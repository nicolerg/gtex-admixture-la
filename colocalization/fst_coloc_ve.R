#!/bin/R

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(lme4)

if(!file.exists('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc_fst-for-regression.RData')){
  
  load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-1e-04-20200414.RData')
  head(master_coloc)
  nrow(unique(master_coloc[,.(gene_id, tissue)])) # 4529, same as pruned egenes 
  master_coloc[,variant_local := sapply(lead_variants_local, function(x) unname(unlist(strsplit(x, ';')))[1])]
  master_coloc[,variant_global := sapply(lead_variants_global, function(x) unname(unlist(strsplit(x, ';')))[1])]
  master_coloc[,pos_global := sapply(variant_global, function(x) as.numeric(unname(unlist(strsplit(x, '_')))[2]))]
  master_coloc[,pos_local := sapply(variant_local, function(x) as.numeric(unname(unlist(strsplit(x, '_')))[2]))]
  
  master_coloc[,chr_global := sapply(variant_global, function(x) gsub('chr','',unname(unlist(strsplit(x, '_')))[1]))]
  master_coloc[,chr_local := sapply(variant_local, function(x) gsub('chr','',unname(unlist(strsplit(x, '_')))[1]))]
  
  load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/merged-all_no_max.RData')
  merged[,POS := as.numeric(POS)]
  
  # subset merged down to the variants we need to make things easier
  gsub = master_coloc[,.(chr_global, pos_global)]
  lsub = master_coloc[,.(chr_local, pos_local)]
  colnames(gsub) = colnames(lsub) = c('CHROM','POS')
  vsub = rbindlist(list(gsub,lsub))
  vsub = unique(vsub)
  vsub = merge(vsub, merged, by=c('CHROM','POS'))
  
  head(vsub)
  # convert to numeric
  mycols = colnames(vsub)[grepl('AFR_|EUR_',colnames(vsub))]
  vsub[,(mycols) := lapply(.SD, as.numeric), .SDcols = mycols]
  
  # get max within each continent
  mycols = colnames(vsub)[grepl('AFR_',colnames(vsub))]
  vsub[,max_fst_AFR := max(abs(.SD),na.rm=TRUE), by=c('CHROM','POS'), .SDcols = mycols]
  vsub[,which_max_fst_AFR := colnames(.SD)[which.max(abs(.SD))], by=c('CHROM','POS'), .SDcols = mycols]
  mycols = colnames(vsub)[grepl('EUR_',colnames(vsub))]
  vsub[,max_fst_EUR := max(abs(.SD),na.rm=TRUE), by=c('CHROM','POS'), .SDcols = mycols]
  vsub[,which_max_fst_EUR := colnames(.SD)[which.max(abs(.SD))], by=c('CHROM','POS'), .SDcols = mycols]
  
  # get max between continent fst 
  afr_eur_fst = fread('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/EUR_AFR.weir.fst',sep='\t',header=TRUE)
  afr_eur_fst = afr_eur_fst[WEIR_AND_COCKERHAM_FST == '-nan', WEIR_AND_COCKERHAM_FST := 0]
  afr_eur_fst = afr_eur_fst[WEIR_AND_COCKERHAM_FST != 'WEIR_AND_COCKERHAM_FST']
  afr_eur_fst[,WEIR_AND_COCKERHAM_FST := as.numeric(WEIR_AND_COCKERHAM_FST)]
  afr_eur_fst[,POS := as.numeric(POS)]
  
  vsub = vsub[,.(CHROM,POS,max_fst_AFR,which_max_fst_AFR,max_fst_EUR,which_max_fst_EUR)]
  vmerged = merge(vsub, afr_eur_fst, by=c('CHROM','POS'), all.x=T)
  setnames(vmerged, "WEIR_AND_COCKERHAM_FST", "AFR_EUR_fst")
  vmerged[max_fst_EUR == -Inf, max_fst_EUR := NA]
  vmerged[max_fst_AFR == -Inf, max_fst_AFR := NA]
  
  # get max within 
  vmerged[,max_within_fst := max(abs(.SD),na.rm=TRUE), by=c('CHROM','POS'), .SDcols = c('max_fst_AFR','max_fst_EUR')]
  
  # merge back with egenes 
  
  masterg = merge(master_coloc, vmerged, by.x=c('chr_global','pos_global'), by.y=c('CHROM','POS'), all.x = T)
  masterl = merge(master_coloc, vmerged, by.x=c('chr_local','pos_local'), by.y=c('CHROM','POS'), all.x = T)
  master_coloc_fst = merge(masterg, masterl, by=colnames(master_coloc), suffixes = c('_global','_local'))
  
  master_coloc_fst[,c('variant_local','variant_global','pos_global','pos_local','chr_global','chr_local') := NULL]
  
  save(master_coloc_fst, file='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc_fst-for-regression.RData')
}else{
  load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc_fst-for-regression.RData')
}

# try some different regressions 

# first, the original one in the manuscript 
# get max coloc per gene per tissue per method
# m1_coloc = master_coloc_fst[,list(coloc_h4_global = max(coloc_h4_global),
#                             coloc_h4_local = max(coloc_h4_local),
#                             gwas_trait_local = gwas_trait[which.max(coloc_h4_local)],
#                             gwas_trait_global = gwas_trait[which.max(coloc_h4_global)]), by=.(gene_id, tissue, 
#                                                                                   max_within_fst_global, 
#                                                                                   max_within_fst_local, 
#                                                                                   max_fst_AFR_global,
#                                                                                   max_fst_EUR_global,
#                                                                                   AFR_EUR_fst_global,
#                                                                                   max_fst_AFR_local,
#                                                                                   max_fst_EUR_local,
#                                                                                   AFR_EUR_fst_local,
#                                                                                   pval_nominal_global, 
#                                                                                   pval_nominal_local)]
# head(m1_coloc)
# 
# summary(lm(coloc_h4_local ~ max_within_fst_local + AFR_EUR_fst_local, data=m1_coloc))
# summary(lm(coloc_h4_global ~ max_within_fst_global + AFR_EUR_fst_global, data=m1_coloc))

# now, with a linear mixed model and *all* colocalization tests 
head(master_coloc_fst)

# rearrange into narrower df
local = master_coloc_fst[,.(gene_id, tissue, pval_nominal_local, lead_variants_local, ref_snp, gwas_trait, n_snps_coloc, coloc_h4_local, 
                            n_snps_finemap, finemap_clpp_local, `-log_gwas_pval_local`, `-log_eqtl_pval_local`, 
                            max_fst_AFR_local, which_max_fst_AFR_local, max_fst_EUR_local, which_max_fst_EUR_local, AFR_EUR_fst_local, max_within_fst_local)]
local[,method := 'local']

global = master_coloc_fst[,.(gene_id, tissue, pval_nominal_global, lead_variants_global, ref_snp, gwas_trait, n_snps_coloc, coloc_h4_global, 
                            n_snps_finemap, finemap_clpp_global, `-log_gwas_pval_global`, `-log_eqtl_pval_global`, 
                            max_fst_AFR_global, which_max_fst_AFR_global, max_fst_EUR_global, which_max_fst_EUR_global, AFR_EUR_fst_global, max_within_fst_global)]
global[,method := 'global']

colnames(local) = colnames(global) = c('gene_id','tissue','pval_nominal','lead_variants','ref_snp','gwas_trait','n_snps_coloc','coloc_h4',
                                       'n_snps_finemap','finemap_clpp','log10_gwas_pval','log10_eqtl_pval','max_fst_AFR','which_max_fst_AFR',
                                       'max_fst_EUR','which_max_fst_EUR','AFR_EUR_fst','max_within_fst','method')

master_coloc_fst_narrow = rbindlist(list(local,global))
head(master_coloc_fst_narrow)

# lmm <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method*AFR_EUR_fst + method*max_within_fst,
#             data = master_coloc_fst_narrow)
# summary(lmm)

# what if I only include good colocs?
subset = master_coloc_fst[coloc_h4_global > 0.5 | coloc_h4_local > 0.5, .(gene_id, tissue, gwas_trait)]
m2_coloc = merge(master_coloc_fst_narrow, subset, by=c('gene_id','tissue','gwas_trait'))
# lmm <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst,
#             data = m2_coloc)
# summary(lmm)

# # are between_fst and within_fst correlated?
# cor(m2_coloc[,AFR_EUR_fst], m2_coloc[,max_within_fst], use="complete.obs") # 0.14
# cor(m2_coloc[,AFR_EUR_fst], m2_coloc[,max_fst_EUR], use="complete.obs") # 0.226
# cor(m2_coloc[,AFR_EUR_fst], m2_coloc[,max_fst_AFR], use="complete.obs") # -0.041

# TAH DAHHHHHHH 
####################################
m3_coloc = m2_coloc[,.(tissue, gene_id, gwas_trait, coloc_h4, AFR_EUR_fst, max_within_fst, method)]
m3_coloc = m3_coloc[complete.cases(m3_coloc)]
lmm <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
            data = m3_coloc)
summary(lmm)
# Linear mixed model fit by REML ['lmerMod']
# Formula: coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) +  
#   method + AFR_EUR_fst + max_within_fst + method:max_within_fst +      method:AFR_EUR_fst
# Data: m3_coloc
# 
# REML criterion at convergence: -1891.4
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.3963 -0.3957  0.0256  0.4078  4.1601 
# 
# Random effects:
#   Groups     Name        Variance  Std.Dev.
# gene_id    (Intercept) 0.0116678 0.10802 
# gwas_trait (Intercept) 0.0034346 0.05861 
# tissue     (Intercept) 0.0007913 0.02813 
# Residual               0.0071114 0.08433 
# Number of obs: 1155, groups:  gene_id, 207; gwas_trait, 89; tissue, 7
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                 0.59459    0.01776  33.476
# methodlocal                -0.04847    0.01130  -4.290
# AFR_EUR_fst                 0.06767    0.03448   1.963
# max_within_fst             -0.95998    0.17342  -5.536
# methodlocal:max_within_fst  1.09245    0.22018   4.962
# methodlocal:AFR_EUR_fst     0.01590    0.02858   0.556
# 
# Correlation of Fixed Effects:
#   (Intr) mthdlc AFR_EU mx_wt_ mth:__
# methodlocal -0.226                            
# AFR_EUR_fst -0.287  0.028                     
# mx_wthn_fst -0.343  0.370 -0.137              
# mthdlcl:m__  0.180 -0.758  0.195 -0.618       
# mt:AFR_EUR_  0.081 -0.317 -0.364  0.132 -0.177
####################################

# method:max_within_fst
lmm0 <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:AFR_EUR_fst,
            data = m3_coloc)
anova(lmm0,lmm,refit=T)
# p-value: 8.373e-07

# methodlocal
lmm0 <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + AFR_EUR_fst + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
            data = m3_coloc)
anova(lmm0,lmm,refit=T)
# p-value: 1.947e-05

# AFR_EUR_fst
lmm0 <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
            data = m3_coloc)
anova(lmm0,lmm,refit=T)
# p-value < 2.2e-16

# max_within_fst
lmm0 <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + method:max_within_fst + method:AFR_EUR_fst,
            data = m3_coloc)
anova(lmm0,lmm,refit=T)
# p-value < 2.2e-16

# methodlocal:AFR_EUR_fst
lmm0 <- lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:max_within_fst,
            data = m3_coloc)
anova(lmm0,lmm,refit=T)
# p-value: 0.576

####################################

# plot rediduals 
resids = data.table(resid=resid(lmm), AFR_EUR_fst=m3_coloc[,AFR_EUR_fst], max_within_fst=m3_coloc[,max_within_fst], method=m3_coloc[,method])
head(resids)
resids[method=='global',method := 'GlobalAA']
resids[method=='local',method := 'LocalAA']
g1 = ggplot(resids, aes(x=AFR_EUR_fst, y=resid)) + 
  geom_point() + 
  facet_wrap(~method) +
  theme_classic() +
  labs(x=bquote('F'['ST'*","*'between']), y='LMM residual')
g2 = ggplot(resids, aes(x=max_within_fst, y=resid)) + 
  geom_point() + 
  facet_wrap(~method) +
  theme_classic() +
  labs(x=bquote('F'['ST'*","*'within']), y='LMM residual') 

pdf('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/plots/lmm_resid.pdf', width=5.5,height=5)
grid.arrange(g1,g2,
             layout_matrix=rbind(c(1),
                                 c(2))
)
dev.off()

# does finemap give me the same result?
# first try with best finemap coloc
subset = unique(master_coloc_fst[!is.na(finemap_clpp_local) & !is.na(finemap_clpp_global), .(gene_id, tissue, gwas_trait)])
m4_finemap = merge(master_coloc_fst_narrow, subset, by=c('gene_id','tissue','gwas_trait'))
m4_finemap[,finemap_clpp := -log10(finemap_clpp)]
m5_finemap = m4_finemap[,.(tissue, gene_id, gwas_trait, finemap_clpp, AFR_EUR_fst, max_within_fst, method)]
m5_finemap = m5_finemap[complete.cases(m5_finemap)]
m5_finemap = m5_finemap[finemap_clpp != Inf]
lmm <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method*AFR_EUR_fst + method*max_within_fst,
            data = m5_finemap)
summary(lmm)
# no, don't get the same result 

# Linear mixed model fit by REML ['lmerMod']
# Formula: finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) +      method * AFR_EUR_fst + method * max_within_fst
# Data: m5_finemap
# 
# REML criterion at convergence: -391201.7
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -16.208  -0.276   0.021   0.353  46.540 
# 
# Random effects:
#   Groups     Name        Variance  Std.Dev.
# gene_id    (Intercept) 1.970e-02 0.140360
# gwas_trait (Intercept) 1.584e-03 0.039805
# tissue     (Intercept) 2.075e-05 0.004555
# Residual               6.267e-03 0.079166
# Number of obs: 183277, groups:  gene_id, 3751; gwas_trait, 142; tissue, 7
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                 3.3571571  0.0044782 749.663
# methodlocal                 0.0004882  0.0008925   0.547
# AFR_EUR_fst                -0.0092623  0.0024438  -3.790
# max_within_fst              0.0239220  0.0136182   1.757
# methodlocal:AFR_EUR_fst     0.0051127  0.0029145   1.754 NS 
# methodlocal:max_within_fst -0.0310976  0.0172025  -1.808 NS 
# 
# Correlation of Fixed Effects:
#   (Intr) mthdlc AFR_EU mx_wt_ m:AFR_
# methodlocal -0.100                            
# AFR_EUR_fst -0.089  0.314                     
# mx_wthn_fst -0.113  0.476 -0.004              
# mt:AFR_EUR_  0.050 -0.532 -0.558 -0.004       
# mthdlcl:m__  0.075 -0.736 -0.009 -0.656  0.004

# calculate p-values 

# method:max_within_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:AFR_EUR_fst,
             data = m5_finemap)
anova(lmm0,lmm,refit=T)
# p-value: 0.07059

# methodlocal
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + AFR_EUR_fst + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
             data = m5_finemap)
anova(lmm0,lmm,refit=T)
# p-value: can't compute

# AFR_EUR_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
             data = m5_finemap)
anova(lmm0,lmm,refit=T)
# p-value: can't compute

# max_within_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + method:max_within_fst + method:AFR_EUR_fst,
             data = m5_finemap)
anova(lmm0,lmm,refit=T)
# p-value: < 2.2e-16

# methodlocal:AFR_EUR_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:max_within_fst,
             data = m5_finemap)
anova(lmm0,lmm,refit=T)
# p-value: 0.07942


# how many signif COLOC overlap with finemap?
m3_coloc[,hit_name := paste(tissue, gene_id, gwas_trait, sep=':')]
master_coloc_fst_narrow[,hit_name := paste(tissue, gene_id, gwas_trait, sep=':')]
coloc_hit_subset = master_coloc_fst_narrow[hit_name %in% m3_coloc[,hit_name]]
m6_finemap = coloc_hit_subset[,.(tissue, gene_id, gwas_trait, finemap_clpp, AFR_EUR_fst, max_within_fst, method)]
m6_finemap = m6_finemap[complete.cases(m6_finemap)]
m6_finemap = m6_finemap[finemap_clpp != Inf]
m6_finemap[,finemap_clpp := -log10(finemap_clpp)]
# 846 out of 1155 signif coloc hits were tested with finemap 
lmm <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method*AFR_EUR_fst + method*max_within_fst,
            data = m6_finemap)
summary(lmm)
# Linear mixed model fit by REML ['lmerMod']
# Formula: finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) +      method * AFR_EUR_fst + method * max_within_fst
# Data: m6_finemap
# 
# REML criterion at convergence: -795.4
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -6.1622 -0.2702  0.0562  0.2943  5.4355 
# 
# Random effects:
#   Groups     Name        Variance  Std.Dev.
# gene_id    (Intercept) 0.0599664 0.24488 
# gwas_trait (Intercept) 0.0166950 0.12921 
# tissue     (Intercept) 0.0005107 0.02260 
# Residual               0.0091954 0.09589 
# Number of obs: 846, groups:  gene_id, 192; gwas_trait, 87; tissue, 7
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                 3.06720    0.02861 107.220
# methodlocal                 0.01377    0.01559   0.883
# AFR_EUR_fst                -0.06455    0.04860  -1.328
# max_within_fst              0.53521    0.23053   2.322
# methodlocal:AFR_EUR_fst     0.03432    0.03833   0.895
# methodlocal:max_within_fst -0.58637    0.28873  -2.031
# 
# Correlation of Fixed Effects:
#   (Intr) mthdlc AFR_EU mx_wt_ m:AFR_
# methodlocal -0.169                            
# AFR_EUR_fst -0.272  0.008                     
# mx_wthn_fst -0.290  0.331 -0.095              
# mt:AFR_EUR_  0.073 -0.385 -0.329  0.099       
# mthdlcl:m__  0.138 -0.767  0.190 -0.581 -0.095

# no, still don't get the same result 

# calculate p-values 

# method:max_within_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:AFR_EUR_fst,
             data = m6_finemap)
anova(lmm0,lmm,refit=T)
# p-value: 0.04155

# methodlocal
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + AFR_EUR_fst + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
             data = m6_finemap)
anova(lmm0,lmm,refit=T)
# p-value: 0.3742

# AFR_EUR_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
             data = m6_finemap)
anova(lmm0,lmm,refit=T)
# p-value: 1

# max_within_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + method:max_within_fst + method:AFR_EUR_fst,
             data = m6_finemap)
anova(lmm0,lmm,refit=T)
# p-value: 1

# methodlocal:AFR_EUR_fst
lmm0 <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:max_within_fst,
             data = m6_finemap)
anova(lmm0,lmm,refit=T)
# p-value: 0.3684


