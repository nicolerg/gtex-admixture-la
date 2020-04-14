library(data.table)
library(fitdistrplus)

timestamp = function(message){
	writeLines(sprintf('%s	%s', as.character(Sys.time()), message))
}

args = commandArgs(trailingOnly=T)
allp = args[1] # Adipose_Subcutaneous-local-global-permuted-chr22-14.tsv
allp = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/Adipose_Subcutaneous/chr22/Adipose_Subcutaneous-local-global-permuted-chr22-14.tsv'

egenes_master = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/egenes_master.RData'
load(egenes_master)

Tissue = unname(unlist(strsplit(basename(allp), '-')))[1]

gc(verbose=F)

timestamp('Find minimum p-value per permutation')

# read from outfile; saved lists didn't work 
permuted_pvalues = fread(allp, sep='\t', header=T)
minimum_pvalues = permuted_pvalues[,list(global_pval_min=min(global_pval_nominal, na.rm=T),
	local_pval_min=min(local_pval_nominal, na.rm=T),
	local_variant_id=variant_id[which.min(local_pval_nominal)],
	global_variant_id=variant_id[which.min(global_pval_nominal)]),
by=.(tissue, gene_id, permutation)]

timestamp('Fit beta distribution for each gene and approximate adjusted p-value')
dt = data.table(gene_id = unique(minimum_pvalues[,gene_id]), 
				tissue = Tissue,
				pval_nominal = NA_real_, 
				shape1=NA_real_,
				shape2=NA_real_,
				beta_est_pval = NA_real_)

beta_adjusted = list()
beta_adjusted[['global']] = copy(dt)
beta_adjusted[['LAVA']] = copy(dt)

get_beta_adjusted = function(pval_vector, m, gene){
	ml = fitdist(pval_vector, "beta", method="mle")
	s1 = ml$estimate[['shape1']]
	s2 = ml$estimate[['shape2']]
	# get actual p-value from egenes_master
	true_pval = as.numeric(unname(unlist(egenes_master[gene_id == gene & method == m & tissue == Tissue, pval_nominal])))[1]
	# get adjusted p-value from dbeta
	adjusted_pval = pbeta(true_pval, s1, s2)
	beta_adjusted[[m]][gene_id == gene, pval_nominal := true_pval]
	beta_adjusted[[m]][gene_id == gene, shape1 := s1]
	beta_adjusted[[m]][gene_id == gene, shape2 := s2]
	beta_adjusted[[m]][gene_id == gene, beta_est_pval := adjusted_pval]
	return(beta_adjusted)
}

for (gene in unique(minimum_pvalues[,gene_id])){
	# get p-value vector
	xg = minimum_pvalues[gene_id==gene,global_pval_min]
	xl = minimum_pvalues[gene_id==gene,local_pval_min]
	xl = xl[!is.na(xl)]
	xg = xg[!is.na(xg)]

	beta_adjusted = get_beta_adjusted(xg, 'global', gene)
	beta_adjusted = get_beta_adjusted(xl, 'LAVA', gene)

}

# > head(beta_adjusted[['global']])
#               gene_id               tissue pval_nominal    shape1   shape2
# 1:  ENSG00000213279.2 Adipose_Subcutaneous 5.597297e-04 0.9179512 1687.373
# 2:  ENSG00000235111.1 Adipose_Subcutaneous 1.759154e-06 0.9730726 1694.427
# 3: ENSG00000188511.12 Adipose_Subcutaneous 1.526999e-07 0.7822290 1264.885
# 4:  ENSG00000279345.1 Adipose_Subcutaneous 4.679057e-05 0.9469126 1082.511
# 5:  ENSG00000236867.1 Adipose_Subcutaneous 5.963614e-06 0.9607258 1063.038
# 6: ENSG00000100425.18 Adipose_Subcutaneous 1.622240e-04 1.0375753 1129.069
#    beta_est_pval
# 1:   0.647743013
# 2:   0.003519959
# 3:   0.001341744
# 4:   0.059163885
# 5:   0.007834650
# 6:   0.154231802

beta_adjusted[['global']][,method := 'global']
beta_adjusted[['LAVA']][,method := 'LAVA']
all_beta_adj = rbindlist(beta_adjusted)

outfile = sprintf('%s.beta_est.tsv',gsub('\\..*', '',allp))
write.table(all_beta_adj, outfile, sep='\t', col.names=T, row.names=F, quote=F)

q()







# practice plotting
tissuecols <- c(Adipose_Subcutaneous="#FF6600",
	Adipose_Visceral_Omentum="#FFAA00",
	Adrenal_Gland="#33DD33",
	Artery_Aorta="#FF5555",
	Artery_Coronary="#FFAA99",
	Artery_Tibial="#FF0000",
	Brain_Amygdala="#EEEE00",
	Brain_Anterior_cingulate_cortex_BA24="#EEEE00",
	Brain_Caudate_basal_ganglia="#EEEE00",
	Brain_Cerebellar_Hemisphere="#EEEE00",
	Brain_Cerebellum="#EEEE00",
	Brain_Cortex="#EEEE00",
	Brain_Frontal_Cortex_BA9="#EEEE00",
	Brain_Hippocampus="#EEEE00",
	Brain_Hypothalamus="#EEEE00",
	Brain_Nucleus_accumbens_basal_ganglia="#EEEE00",
	Brain_Putamen_basal_ganglia="#EEEE00",
	"Brain_Spinal_cord_cervical_c-1"="#EEEE00",
	Brain_Substantia_nigra="#EEEE00",
	Breast_Mammary_Tissue="#33CCCC",
	"Cells_EBV-transformed_lymphocytes"="#CC66FF",
	Cells_Cultured_fibroblasts="#AAEEFF",
	Colon_Sigmoid="#EEBB77",
	Colon_Transverse="#CC9955",
	Esophagus_Gastroesophageal_Junction="#8B7355",
	Esophagus_Mucosa="#552200",
	Esophagus_Muscularis="#BB9988",
	Heart_Atrial_Appendage="#9900FF",
	Heart_Left_Ventricle="#660099",
	Liver="#AABB66",
	Lung="green4", # changed from #9aff00, which is lime green
	Kidney_Cortex="yellow",
	Minor_Salivary_Gland="#99BB88",
	Muscle_Skeletal="#AAAAFF",
	Nerve_Tibial="#FFD700",
	Ovary="#FFAAFF",
	Pancreas="#995522",
	Pituitary="#AAFF99",
	Prostate="#DDDDDD",
	Skin_Not_Sun_Exposed_Suprapubic="#0000FF",
	Skin_Sun_Exposed_Lower_leg="#7777FF",
	Small_Intestine_Terminal_Ileum="#555522",
	Spleen="#778855",
	Stomach="#FFDD99",
	Testis="#AAAAAA",
	Thyroid="#006600",
	Uterus="#FF66FF",
	Vagina="#FF5599",
	Whole_Blood="#FF00BB")

tissuelabs <- c(Adipose_Subcutaneous="Subc. adipose",
	Artery_Tibial="Tibial artery",
	Lung="Lung",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Skin_Not_Sun_Exposed_Suprapubic="NSE skin",
	Whole_Blood='Whole blood')

library(data.table)
library(ggplot2)
all_beta_adj = fread('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/Adipose_Subcutaneous/chr22/Adipose_Subcutaneous-local-global-permuted-chr22-14.beta_est.tsv', sep='\t', header=T)

fits = data.table(method=c('LocalAA','GlobalAA'),
	y=NA_real_,
	slope=NA_real_,
	intercept=NA_real_)

# # linear regression 
# fit_g = lm(I(-log10(beta_est_pval)) ~ I(-log10(pval_nominal)), data=all_beta_adj[method=='global'])
# fit_l = lm(I(-log10(beta_est_pval)) ~ I(-log10(pval_nominal)), data=all_beta_adj[method=='LAVA'])
# 
# fits[method=='LocalAA', slope := fit_l$coefficients[[2]]]
# fits[method=='LocalAA', intercept := fit_l$coefficients[[1]]]
# fits[method=='GlobalAA', slope := fit_g$coefficients[[2]]]
# fits[method=='GlobalAA', intercept := fit_g$coefficients[[1]]]
# fits[,y := intercept + 6*slope]
# 
# all_beta_adj[method == 'LAVA', method:='LocalAA']
# all_beta_adj[method == 'global', method:='GlobalAA']
# 
# ggplot(all_beta_adj, aes(x=-log10(pval_nominal), y=-log10(beta_est_pval), colour=tissue)) +
# 	geom_vline(xintercept=6, linetype='dashed', colour='gray') +
# 	geom_hline(data=fits, aes(yintercept=y), linetype='dashed', colour='gray') +
# 	geom_abline(data=fits, aes(slope=slope, intercept=intercept), linetype='dashed', colour='black') + 
# 	geom_point() + 
# 	theme_classic() + 
# 	labs(x='Nominal p-value (-log10)', y='Beta-estimated adjusted p-value (-log10)') +
# 	facet_wrap(~method) +
#   scale_colour_manual(name='Tissue',values=tissuecols,labels=tissuelabs)


# I don't want methods to be in separate regressions
fit = lm(I(-log10(beta_est_pval)) ~ I(-log10(pval_nominal)), data=all_beta_adj)
slope = fit$coefficients[[2]]
intercept = fit$coefficients[[1]]
y = intercept + 6*slope

all_beta_adj[method == 'LAVA', method:='LocalAA']
all_beta_adj[method == 'global', method:='GlobalAA']
max_beta = ceiling(-log10(min(all_beta_adj[,beta_est_pval])))

g1 = ggplot(all_beta_adj, aes(x=-log10(pval_nominal), y=-log10(beta_est_pval), colour=tissue, shape=method)) +
  geom_vline(xintercept=6, linetype='dashed', colour='gray') +
  geom_hline(yintercept=y, linetype='dashed', colour='gray') +
  geom_abline(slope=slope, intercept=intercept, linetype='dashed', colour='black') + 
  geom_point() + 
  theme_classic() + 
  labs(x='Nominal P-value (-log10)', y='Beta-estimated adjusted P-value (-log10)') +
  scale_colour_manual(name='Tissue',values=tissuecols,labels=tissuelabs) +
  scale_shape_discrete(name='Method') +
  scale_y_continuous(breaks=c(seq(from=0,to=max_beta,by=1), y),
                     labels=as.character(c(seq(from=0,to=max_beta,by=1), round(y,2)))) +
  theme(legend.position = 'none') +
  annotate("text",x=min(-log10(all_beta_adj[,pval_nominal])), y=y+0.02, label=sprintf('adj. P-value = %s', signif(10^-y, 2)),
           vjust=0, hjust=0)

library(qvalue)
global_qval = qvalue(all_beta_adj[method=='GlobalAA',beta_est_pval],pi0 = 1)$qvalues
local_qval = qvalue(all_beta_adj[method=='LocalAA',beta_est_pval],pi0 = 1)$qvalues
all_beta_adj[method=='GlobalAA',qvalue := global_qval]
all_beta_adj[method=='LocalAA',qvalue := local_qval]

# get qvalue corresponding with 1e-6
fit = lm(I(-log10(qvalue)) ~ I(-log10(pval_nominal)), data=all_beta_adj)
slope = fit$coefficients[[2]]
intercept = fit$coefficients[[1]]
y = intercept + 6*slope

max_qval = ceiling(-log10(min(all_beta_adj[,qvalue])))

g2 = ggplot(all_beta_adj, aes(x=-log10(pval_nominal), y=-log10(qvalue), colour=tissue, shape=method)) +
  geom_vline(xintercept=6, linetype='dashed', colour='gray') +
  geom_hline(yintercept=y, linetype='dashed', colour='gray') +
  geom_abline(slope=slope, intercept=intercept, linetype='dashed', colour='black') + 
  geom_point() + 
  theme_classic() + 
  labs(x='Nominal P-value (-log10)', y='qvalue (-log10)') +
  scale_colour_manual(name='Tissue',values=tissuecols,labels=tissuelabs) +
  scale_shape_discrete(name='Method') +
  scale_y_continuous(breaks=c(seq(from=0,to=max_beta,by=0.5), y),
                     labels=as.character(round(c(seq(from=0,to=max_beta,by=0.5), y), 2))) +
  annotate("text",x=min(-log10(all_beta_adj[,pval_nominal])), y=y+0.02, label=sprintf('qvalue = %s', signif(10^-y, 2)),
           vjust=0, hjust=0)

# print both plots

library(grid)
library(gridExtra)

pdf('~/test-adjusted-pvals.pdf',width=8,height=3.5)
grid.arrange(g1, g2,
             layout_matrix=rbind(c(1,1,2,2,2),
                                 c(1,1,2,2,2))
)
dev.off()
