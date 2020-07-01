library(data.table)
library(fitdistrplus)

timestamp = function(message){
	writeLines(sprintf('%s	%s', as.character(Sys.time()), message))
}

indir = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/chunks'

# load egenes_master
load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/egenes_master.RData')

min_all = list()
i = 1
for (chunk in list.files(path=indir, full.names=T)){
  timestamp(chunk)
  permuted_pvalues = fread(chunk, sep='\t', header=F)
  colnames(permuted_pvalues) = c('tissue','gene_id','permutation','variant_id','global_pval_nominal','local_pval_nominal')
  timestamp('Find minimum p-value per permutation')
  minimum_pvalues = permuted_pvalues[,list(global_pval_min=min(global_pval_nominal, na.rm=T),
					   local_pval_min=min(local_pval_nominal, na.rm=T),
					   local_variant_id=variant_id[which.min(local_pval_nominal)],
					   global_variant_id=variant_id[which.min(global_pval_nominal)]),
				     by=.(tissue, gene_id, permutation)]
	min_all[[i]] = minimum_pvalues
	i = i + 1
}

minimum_pvalues = rbindlist(min_all)
head(minimum_pvalues)

save(minimum_pvalues, file='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/minimum_pvalues.RData')

timestamp('Fit beta distribution for each gene and approximate adjusted p-value')
dt = unique(minimum_pvalues[,.(tissue, gene_id)])
dt[,pval_nominal := NA_real_]
dt[,shape1:=NA_real_]
dt[,shape2:=NA_real_]
dt[,beta_est_pval := NA_real_]

beta_adjusted = list()
beta_adjusted[['global']] = copy(dt)
beta_adjusted[['LAVA']] = copy(dt)

get_beta_adjusted = function(pval_vector, m, gene, Tissue){
	ml = fitdist(pval_vector, "beta", method="mle")
	s1 = ml$estimate[['shape1']]
	s2 = ml$estimate[['shape2']]
	# get actual p-value from egenes_master
	true_pval = as.numeric(unname(unlist(egenes_master[gene_id == gene & method == m & tissue == Tissue, pval_nominal])))[1]
	# get adjusted p-value from dbeta
	adjusted_pval = pbeta(true_pval, s1, s2)
	beta_adjusted[[m]][gene_id == gene & tissue == Tissue, pval_nominal := true_pval]
	beta_adjusted[[m]][gene_id == gene & tissue == Tissue, shape1 := s1]
	beta_adjusted[[m]][gene_id == gene & tissue == Tissue, shape2 := s2]
	beta_adjusted[[m]][gene_id == gene & tissue == Tissue, beta_est_pval := adjusted_pval]
	return(beta_adjusted)
}

t = unique(minimum_pvalues[,.(tissue, gene_id)])
for (i in 1:nrow(t)){
  Tissue = t[i, tissue]
  gene = t[i, gene_id]
	# get p-value vector
	xg = minimum_pvalues[gene_id==gene & tissue==Tissue,global_pval_min]
	xl = minimum_pvalues[gene_id==gene & tissue==Tissue,local_pval_min]
	xl = xl[!is.na(xl)]
	xg = xg[!is.na(xg)]

	beta_adjusted = get_beta_adjusted(xg, 'global', gene, Tissue)
	beta_adjusted = get_beta_adjusted(xl, 'LAVA', gene, Tissue)

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

outfile = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/chr22_beta_adjusted_pvalues.tsv'
write.table(all_beta_adj, outfile, sep='\t', col.names=T, row.names=F, quote=F)
q()

# plot
source('/oak/stanford/groups/smontgom/nicolerg/src/gtex-admixture-la/plotting/august_figures.R')

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(qvalue)

all_beta_adj = fread('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/chr22_beta_adjusted_pvalues.tsv', sep='\t', header=T)

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
  scale_y_continuous(breaks=c(seq(from=0,to=max_beta,by=2), y),
                     labels=as.character(c(seq(from=0,to=max_beta,by=2), round(y,2)))) +
  scale_x_continuous(breaks=c(seq(5,15,5),6)) +
  theme(legend.position = 'none') +
  annotate("text",x=7, y=y+0.2, label=sprintf('adj. P-value = %s', signif(10^-y, 2)),
           vjust=0, hjust=0)


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
  scale_y_continuous(breaks=c(seq(from=0,to=max_beta,by=1), y),
                     labels=as.character(round(c(seq(from=0,to=max_beta,by=1), y), 2))) +
  scale_x_continuous(breaks=c(seq(5,15,5),6)) +
  annotate("text",x=9, y=y+0.2, label=sprintf('qvalue = %s', signif(10^-y, 2)),
           vjust=0, hjust=0)

# print both plots
pdf('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/plots/adjusted_p_values.pdf',width=8,height=3.5)
grid.arrange(g1, g2,
             layout_matrix=rbind(c(1,1,2,2,2),
                                 c(1,1,2,2,2))
)
dev.off()
