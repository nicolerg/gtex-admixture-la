#!/bin/R

# clean up coloc results 

library(data.table)

indir = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC'
setwd(indir)

fixed = fread(cmd=sprintf('zcat cc_fixed_coloc.tsv.gz'), sep='\t', header=T)
nrow(fixed)
head(fixed)
# make tissue column
fixed[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
fixed[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
table(fixed[,tissue])
table(fixed[,method])
fixed[,eqtl_file := NULL]
fixed[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
fixed[,c('base_gwas_file','gwas_trait'):=NULL]
fixed = fixed[order(tissue, feature, gwas)]

depr = fread(cmd=sprintf('zcat all_coloc_results.tsv.gz'), sep='\t', header=T)
# format 
depr[, tissue := sapply(eqtl_file, function(x) unname(unlist(strsplit(x, '_filtered')))[1])]
depr[, method := sapply(eqtl_file, function(x) ifelse(grepl('global',x), 'global','local'))]
depr[,eqtl_file := NULL]
depr[,gwas := sapply(gwas_trait, function(x) gsub('\\.formatted.*','',x))]
depr[,gwas := sapply(gwas, function(x) gsub('\\.txt\\.gz','',x))]
depr[,gwas := sapply(gwas, function(x) gsub('coloc_imputed_','imputed_',x))]
depr[,c('base_gwas_file','gwas_trait'):=NULL]
depr = depr[order(tissue, feature, gwas)]

# now remove gwas present in "fixed"
gwas_fixed = unique(fixed[,gwas])
length(unique(depr[,gwas]))
depr = depr[!gwas %in% gwas_fixed]
merged = rbind(depr, fixed)

merged = merged[order(tissue, feature, gwas)]
table(merged[,gwas])
table(merged[,method], merged[,tissue])

head(merged)
setnames(merged, "gwas", "gwas_trait")

# split local and global 
l = merged[method=='local']
g = merged[method=='global']
l = unique(l)
g = unique(g)

l[,method := NULL]
g[,method := NULL]

merged = merge(g, l, by=c('ref_snp','n_snps','tissue','gwas_trait','feature'), suffixes=c('_global','_local'), all=T)
write.table(merged, 'merged_coloc_results.tsv', sep='\t', col.names=T, row.names=F)

# how many genes are excluded if we remove incomplete cases?
incomplete = merged[!complete.cases(merged)]
incomplete = incomplete[,.(tissue, feature, gwas_trait)]
incomplete = unique(incomplete)
nrow(incomplete)
merged_complete = merged[complete.cases(merged)]
rows = c()
for (i in 1:nrow(incomplete)){
  if(nrow(merged_complete[tissue == incomplete[i,tissue] & feature == incomplete[i,feature] & gwas_trait == incomplete[i,gwas_trait]]) == 0){
    rows = c(rows, i)
  }
}
incomplete = incomplete[rows]
# are these genes I even wanted to test?
inc = unique(incomplete[,.(tissue, feature)])

load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/egenes_master-20200326.RData')
egenes = egenes_master[pval_nominal_global < 1e-4 & pval_nominal_local < 1e-4]
egenes = egenes[overlapping_lead_variants==0]

missing = 0
for (i in 1:nrow(inc)){
  if(nrow(egenes[tissue == inc[i,tissue] & gene_id == inc[i,feature]]) > 0){
    missing = missing + 1
  }
}
missing # I have no idea why these loci were skipped for only one adjustment method (N = 38)

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

load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/egenes_master-20200326.RData')
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

save(master_coloc, file='master_coloc-1e-04-20200410.RData')
