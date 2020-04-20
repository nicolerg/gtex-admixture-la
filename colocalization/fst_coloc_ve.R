#!/bin/R

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(lmerTest)

if(!file.exists('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc_fst-for-regression_20200418.RData')){
  
  load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-1e-04-20200418.RData')
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
  
  save(master_coloc_fst, file='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc_fst-for-regression_20200418.RData')
}else{
  load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc_fst-for-regression_20200418.RData')
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
local = master_coloc_fst[,.(gene_id, tissue, pval_nominal_local, lead_variants_local, gwas_trait, ref_snp_local_coloc, n_snps_local_coloc, coloc_h4_local, 
                            ref_snp_local_finemap, n_snps_global_finemap, finemap_clpp_local, gwas_pval_log10_local, eqtl_pval_log10_local, 
                            max_fst_AFR_local, which_max_fst_AFR_local, max_fst_EUR_local, which_max_fst_EUR_local, AFR_EUR_fst_local, max_within_fst_local)]
local[,method := 'local']

global = master_coloc_fst[,.(gene_id, tissue, pval_nominal_global, lead_variants_global, gwas_trait, ref_snp_global_coloc, n_snps_global_coloc, coloc_h4_global, 
                             ref_snp_global_finemap, n_snps_global_finemap, finemap_clpp_global, gwas_pval_log10_global, eqtl_pval_log10_global,
                            max_fst_AFR_global, which_max_fst_AFR_global, max_fst_EUR_global, which_max_fst_EUR_global, AFR_EUR_fst_global, max_within_fst_global)]
global[,method := 'global']

colnames(local) = colnames(global) = c('gene_id','tissue','pval_nominal','lead_variants','gwas_trait','ref_snps_coloc','n_snps_coloc','coloc_h4',
                                       'ref_snp_finemap','n_snps_finemap','finemap_clpp','gwas_pval_log10_local','eqtl_pval_log10_local','max_fst_AFR','which_max_fst_AFR',
                                       'max_fst_EUR','which_max_fst_EUR','AFR_EUR_fst','max_within_fst','method')

master_coloc_fst_narrow = rbindlist(list(local,global))
head(master_coloc_fst_narrow)

m_sub = master_coloc_fst_narrow[,.(max_fst_AFR, max_fst_EUR, AFR_EUR_fst, max_within_fst)]
m_sub = unique(m_sub)

ggplot(m_sub, aes(x=AFR_EUR_fst, y=max_within_fst)) +
  geom_point()

ggplot(m_sub, aes(x=AFR_EUR_fst, y=max_fst_AFR)) +
  geom_point()

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
lmm <- lmerTest::lmer(coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:max_within_fst + method:AFR_EUR_fst,
            data = m3_coloc)
summary(lmm)
coloc = data.table(summary(lmm)$coefficients, keep.rownames = T)
colnames(coloc) = c('coefficient','estimate','se','df','t','pvalue')
coloc = coloc[coefficient != '(Intercept)']

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: coloc_h4 ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) +  
#     method + AFR_EUR_fst + max_within_fst + method:max_within_fst +      method:AFR_EUR_fst
#    Data: m3_coloc
# 
# REML criterion at convergence: -2467
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.6930 -0.4362  0.0375  0.4663  3.8519 
# 
# Random effects:
#  Groups     Name        Variance Std.Dev.
#  gene_id    (Intercept) 0.010894 0.10438 
#  gwas_trait (Intercept) 0.004946 0.07033 
#  tissue     (Intercept) 0.001559 0.03948 
#  Residual               0.006738 0.08209 
# Number of obs: 1438, groups:  gene_id, 231; gwas_trait, 91; tissue, 7
# 
# Fixed effects:
#                              Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                 5.955e-01  2.028e-02  1.414e+01  29.368 4.43e-14 ***
# methodlocal                -2.070e-02  9.918e-03  1.115e+03  -2.087 0.037091 *  
# AFR_EUR_fst                -9.879e-03  2.862e-02  1.230e+03  -0.345 0.729982    
# max_within_fst             -9.345e-01  1.545e-01  1.285e+03  -6.048 1.93e-09 ***
# methodlocal:max_within_fst  6.951e-01  1.923e-01  1.158e+03   3.614 0.000315 ***
# methodlocal:AFR_EUR_fst    -1.030e-02  2.494e-02  1.051e+03  -0.413 0.679582 

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
subset = unique(master_coloc_fst[finemap_clpp_local > 0.01 | finemap_clpp_global > 0.01, .(gene_id, tissue, gwas_trait)])
m4_finemap = merge(master_coloc_fst_narrow, subset, by=c('gene_id','tissue','gwas_trait'))
m5_finemap = m4_finemap[,.(tissue, gene_id, gwas_trait, finemap_clpp, AFR_EUR_fst, max_within_fst, method)]
m5_finemap = m5_finemap[complete.cases(m5_finemap)]
hist(m5_finemap[,finemap_clpp], breaks=500)
lmm <- lmer(finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) + method*AFR_EUR_fst + method*max_within_fst,
            data = m5_finemap)
summary(lmm)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: finemap_clpp ~ (1 | tissue) + (1 | gwas_trait) + (1 | gene_id) +  
#     method * AFR_EUR_fst + method * max_within_fst
#    Data: m5_finemap
# 
# REML criterion at convergence: -1013.5
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.6052 -0.2946 -0.0833  0.1514  6.0145 
# 
# Random effects:
#  Groups     Name        Variance  Std.Dev.
#  gene_id    (Intercept) 0.0032721 0.05720 
#  gwas_trait (Intercept) 0.0000000 0.00000 
#  tissue     (Intercept) 0.0000000 0.00000 
#  Residual               0.0008662 0.02943 
# Number of obs: 306, groups:  gene_id, 106; gwas_trait, 56; tissue, 7
# 
# Fixed effects:
#                              Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                  0.045392   0.008780 261.651039   5.170 4.66e-07 ***
# methodlocal                 -0.015320   0.008048 201.041871  -1.904   0.0584 .  
# AFR_EUR_fst                 -0.029771   0.021791 251.586609  -1.366   0.1731    
# max_within_fst               0.016867   0.138962 237.415139   0.121   0.9035    
# methodlocal:AFR_EUR_fst      0.014994   0.023123 201.368781   0.648   0.5174    
# methodlocal:max_within_fst   0.139331   0.155161 205.156969   0.898   0.3703    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) mthdlc AFR_EU mx_wt_ m:AFR_
# methodlocal -0.403                            
# AFR_EUR_fst -0.406  0.159                     
# mx_wthn_fst -0.615  0.423  0.064              
# mt:AFR_EUR_  0.118 -0.450 -0.380  0.027       
# mthdlcl:m__  0.390 -0.795  0.018 -0.645  0.035
# convergence code: 0
# boundary (singular) fit: see ?isSingular


# remove tissue and GWAS trait so it's not singular 
lmm <- lmer(finemap_clpp ~ (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:AFR_EUR_fst + method:max_within_fst,
            data = m5_finemap)
finemap = data.table(summary(lmm)$coefficients, keep.rownames = T)
colnames(finemap) = c('coefficient','estimate','se','df','t','pvalue')
finemap = finemap[coefficient != '(Intercept)']
anova(lmm)

finemap[,method := 'FINEMAP CLPP']
coloc[,method := 'COLOC PP4']
lmer_coef = rbindlist(list(coloc, finemap))

save(lmer_coef, file='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/lmer_coef-to-plot.RData')

# load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/lmer_coef-to-plot.RData')
# ggplot(lmer_coef, aes(y = coefficient, x = estimate)) +
#   geom_vline(xintercept=0, linetype='dashed') +
#   geom_point() +
#   geom_errorbar(aes(xmin=estimate-se, xmax=estimate+se), width=.2) +
#   theme_classic() +
#   facet_wrap(~method, ncol=2) + 
#   theme(axis.text.y=element_text(colour='black'),
#         axis.title.y=element_blank()) +
#   labs(x='Coefficient estimate') +
#   scale_y_discrete(limits = rev(c('methodlocal',
#                               'AFR_EUR_fst',
#                               'max_within_fst',
#                               'methodlocal:AFR_EUR_fst',
#                               'methodlocal:max_within_fst')),
#                    labels=c('methodlocal:max_within_fst'=expression(LocalAA:F[ST*","*within]),
#                             'methodlocal:AFR_EUR_fst'=expression(LocalAA:F[ST*","*between]),
#                             'methodlocal'='LocalAA',
#                             'max_within_fst'=expression(F[ST*","*within]),
#                             'AFR_EUR_fst'=expression(F[ST*","*between])))


# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: finemap_clpp ~ (1 | gene_id) + method + AFR_EUR_fst + max_within_fst +  
#     method:AFR_EUR_fst + method:max_within_fst
#    Data: m5_finemap
# 
# REML criterion at convergence: -1013.5
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.6052 -0.2946 -0.0833  0.1514  6.0145 
# 
# Random effects:
#  Groups   Name        Variance  Std.Dev.
#  gene_id  (Intercept) 0.0032721 0.05720 
#  Residual             0.0008662 0.02943 
# Number of obs: 306, groups:  gene_id, 106
# 
# Fixed effects:
#                              Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                  0.045392   0.008780 261.651039   5.170 4.66e-07 ***
# methodlocal                 -0.015320   0.008048 201.041871  -1.904   0.0584 .  
# AFR_EUR_fst                 -0.029771   0.021791 251.586609  -1.366   0.1731    
# max_within_fst               0.016867   0.138962 237.415139   0.121   0.9035    
# methodlocal:AFR_EUR_fst      0.014994   0.023123 201.368781   0.648   0.5174    
# methodlocal:max_within_fst   0.139331   0.155161 205.156969   0.898   0.3703    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) mthdlc AFR_EU mx_wt_ m:AFR_
# methodlocal -0.403                            
# AFR_EUR_fst -0.406  0.159                     
# mx_wthn_fst -0.615  0.423  0.064              
# mt:AFR_EUR_  0.118 -0.450 -0.380  0.027       
# mthdlcl:m__  0.390 -0.795  0.018 -0.645  0.035
# > anova(lmm)
# Type III Analysis of Variance Table with Satterthwaite's method
# Sum Sq    Mean Sq NumDF  DenDF F value Pr(>F)  
# method                0.00313852 0.00313852     1 201.04  3.6235 0.0584 .
# AFR_EUR_fst           0.00103061 0.00103061     1 270.04  1.1899 0.2763  
# max_within_fst        0.00056829 0.00056829     1 252.57  0.6561 0.4187  
# method:AFR_EUR_fst    0.00036421 0.00036421     1 201.37  0.4205 0.5174  
# method:max_within_fst 0.00069843 0.00069843     1 205.16  0.8064 0.3703  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# does mixed versus EUR predict anything?
mixed = c('Wojcik',
          'pgc.scz2',
          'RA_OKADA_TRANS_ETHNIC',
          'ISGC_Malik_2016_METASTROKE_all_strokes',
          'ILAE_Genetic_generalised_epilepsy',
          'DIAGRAM_T2D_TRANS_ETHNIC',
          'CARDIoGRAM_C4D_CAD_ADDITIVE')
master_coloc_fst_narrow[,EUR := sapply(gwas_trait, function(x) ifelse(any(sapply(mixed, function(y) grepl(y, x))), 0, 1))]
subset = unique(master_coloc_fst[finemap_clpp_local > 0.01 | finemap_clpp_global > 0.01, .(gene_id, tissue, gwas_trait)])
m4_finemap = merge(master_coloc_fst_narrow, subset, by=c('gene_id','tissue','gwas_trait'))
m5_finemap = m4_finemap[,.(tissue, gene_id, gwas_trait, finemap_clpp, AFR_EUR_fst, max_within_fst, method, EUR)]
m5_finemap = m5_finemap[complete.cases(m5_finemap)]

lmm <- lmer(finemap_clpp ~ (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:AFR_EUR_fst + method:max_within_fst + EUR,
            data = m5_finemap)
summary(lmm)
anova(lmm)

lmm <- lmer(finemap_clpp ~ (1 | gene_id) + method + AFR_EUR_fst + max_within_fst + method:AFR_EUR_fst + method:max_within_fst + EUR*method,
            data = m5_finemap)
summary(lmm)
anova(lmm)


