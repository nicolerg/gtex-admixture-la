library(data.table)

# pairs_for_ld.R and plink_ld.sh

sets = fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/20200325/evariant_pairs_test_ld.txt', sep='\t', header=F)
ld_results = fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/20200325/pairs-out-REVISED.tsv', sep='\t', header=F)

# merge back with sets 
colnames(sets) = c('gene_id','tissue','lead_local','lead_global')
colnames(ld_results) = c('lead_local','lead_global','LD')

ld = merge(sets, ld_results, by=c('lead_local','lead_global'))
ld = unique(ld[,.(gene_id,tissue,LD)])

# merge back with egenes_master
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/merged/egenes_master.RData')

# make egenes_master wide 
egenes_master <- egenes_master[,list(lead_variants=paste(variant_id, collapse=';')),by=c('gene_id','count','pval_nominal','method','tissue')]
#egenes_master = egenes_master[,list(label=paste(gene_id,tissue,sep=':')), by=c('gene_id','count','lead_variants','pval_nominal','method','tissue')]
# 400 duplicated cases when you include slope in "by" -- tied p-values can have opposite signs for effect size
g <- unique(egenes_master[method == 'global'])
l <- unique(egenes_master[method == 'LAVA'])
l[,method:=NULL]
g[,method:=NULL]
egenes_master = merge(g, l, by=c('gene_id','tissue','count'), suffixes=c('_global','_local'))
head(egenes_master)

# add an LD columns
egenes_master = merge(egenes_master, ld, by=c('gene_id','tissue'), all.x=T)
# anything not in "ld" has same lead SNP
egenes_master[, overlapping_lead_variants := ifelse(is.na(LD), 1, 0)]

head(egenes_master)

table(egenes_master[,tissue])
table(egenes_master[,overlapping_lead_variants])

print(nrow(egenes_master[overlapping_lead_variants == 0 & LD <1]))

save(egenes_master, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/merged/egenes_master-20200326.RData')
