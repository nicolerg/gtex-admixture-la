#!/bin/R

# DEPRECATED 

# clean up merged coloc results 

library(data.table)

indir <- 'coloc/output/gtex-admixed-difflead/2019-05-24_14-28-31.840186_117admixed'

my_results <- paste0(indir, '/merged_coloc.txt')
coloc <- fread(my_results, header=TRUE, sep='\t')

# make tissue column
coloc[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]

# clean up trait column
coloc[, gwas_trait := gsub('coloc_imputed_','',gwas_trait)]
coloc[, gwas_trait := gsub('\\.txt.*','',gwas_trait)]
coloc[, base_gwas_file := NULL]

table(coloc[,tissue])

coloc[,method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global', 'local'))]

table(coloc[,method])

# split local and global 
l = coloc[method=='local']
g = coloc[method=='global']

l[,method := NULL]
g[,method := NULL]
l[,eqtl_file := NULL]
g[,eqtl_file := NULL]

merged = merge(g, l, by=c('ref_snp','n_snps','tissue','gwas_trait','feature'), suffixes=c('_global','_local'))

# in some cases, the same gene is tested for coloc twice
# this is an artifact of collaping eGene on both 'pval_nominal' and 'slope' upstream
# cases where two lead p-values were the same but had opposite directions of effect (same magnitude) were tested twice 
# take the coloc where the most SNPs were tested
m = merged[,list(n_snps=max(n_snps),
		ref_snp = ref_snp[which.max(n_snps)],
		clpp_h0_global = clpp_h0_global[which.max(n_snps)],
		clpp_h1_global = clpp_h1_global[which.max(n_snps)],
		clpp_h2_global = clpp_h2_global[which.max(n_snps)],
		clpp_h3_global = clpp_h3_global[which.max(n_snps)],
		clpp_h4_global = clpp_h4_global[which.max(n_snps)],
		clpp_h0_local = clpp_h0_local[which.max(n_snps)],
		clpp_h1_local = clpp_h1_local[which.max(n_snps)],
		clpp_h2_local = clpp_h2_local[which.max(n_snps)],
		clpp_h3_local = clpp_h3_local[which.max(n_snps)],
		clpp_h4_local = clpp_h4_local[which.max(n_snps)]), 
	by = c('tissue','gwas_trait','feature') ]

coloc = m

# use LD (PLINK) to filter out loci where R^2==1

load('merged/egenes_master-20191030.RData')
egenes_master = egenes_master[,.(gene_id,tissue,LD,overlapping_lead_variants,pval_nominal_local,pval_nominal_global)]
egenes_master = egenes_master[pval_nominal_local < 1e-4 & pval_nominal_global < 1e-4]

nrow(coloc)
master_coloc = merge(egenes_master,coloc, by.x=c('gene_id','tissue'), by.y=c('feature','tissue'), all.y=T)
nrow(master_coloc)
master_coloc = master_coloc[LD < 1 & overlapping_lead_variants == 0]
nrow(master_coloc)

head(master_coloc)

# check for duplicates
master_coloc = unique(master_coloc)
nrow(master_coloc)
print(nrow(unique(master_coloc[,.(tissue,gwas_trait,gene_id)])))

save(master_coloc, file='coloc/master_coloc-1e-04-2.RData')

