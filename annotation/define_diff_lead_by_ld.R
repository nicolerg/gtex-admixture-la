library(data.table)

# pairs_for_ld.R and plink_ld.sh

sets = fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/20191030/evariant_pairs_test_ld.txt', sep='\t', header=F)
ld_results = fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/20191030/pairs-out-merged.txt', sep='\t', header=F)

# merge back with sets 
colnames(sets) = c('gene_id','tissue','lead_local','lead_global')
colnames(ld_results) = c('lead_local','lead_global','LD')

ld = merge(sets, ld_results, by=c('lead_local','lead_global'))
ld = unique(ld[,.(gene_id,tissue,LD)])

# merge back with egenes_master
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

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

# # "sets" excludes overlapping gene/tissue pairs 
# sets = sets[,list(include=paste(gene_id,tissue,sep=':')), by=c('gene_id','tissue','lead_local','lead_global')]

# egenes_master[,overlapping_lead_variants := sapply(label, function(x) ifelse(x %in% sets[,include], 0, 1))]
# table(egenes_master[,overlapping_lead_variants])
# egenes_master[,label := NULL]

# add an LD columns
egenes_master = merge(egenes_master, ld, by=c('gene_id','tissue'), all.x=T)
# anything not in "ld" has same lead SNP
egenes_master[, overlapping_lead_variants := ifelse(is.na(LD), 1, 0)]

head(egenes_master)

table(egenes_master[,tissue])
table(egenes_master[,overlapping_lead_variants])

# # compare to original 
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')
# egenes_master = unique(egenes_master[,.(gene_id,tissue,method)])
# table(egenes_master[,tissue], egenes_master[,method])

print(nrow(egenes_master[overlapping_lead_variants == 0 & LD <1]))

save(egenes_master, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master-20191030.RData')

q()












### DEPRECATED

# fix definition of "different lead SNPs"

ld <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/local-global-pairs-out.txt',sep='\t',header=F)
#load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-04.RData')
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/m-r2_by_fst_bin.RData')

colnames(ld) <- c('variant_id_global','variant_id_local','D')

m <- merge(m, ld, by=c('variant_id_global','variant_id_local'),all=T)

nrow(m)
m1 <- m[D < 1]
nrow(m1)
head(m)
writeLines('')

# how is the best way to fix this?

m1 <- m1[,.(gene_id, tissue)]

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-04.RData')
diff_lead_fixed <- merge(diff_lead, m1, by=c('gene_id','tissue'))
nrow(unique(diff_lead_fixed, by=c('gene_id','tissue')))
save(diff_lead_fixed, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_fixed_1e-04.RData')

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-05.RData')
diff_lead_fixed <- merge(diff_lead, m1, by=c('gene_id','tissue'))
nrow(unique(diff_lead_fixed, by=c('gene_id','tissue')))
save(diff_lead_fixed, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_fixed_1e-05.RData')

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-06.RData')
diff_lead_fixed <- merge(diff_lead, m1, by=c('gene_id','tissue'))
nrow(unique(diff_lead_fixed, by=c('gene_id','tissue')))
save(diff_lead_fixed, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_fixed_1e-06.RData')

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')
head(egenes_master)

m2 <- m[,.(gene_id, tissue, D)]
head(m2)

egenes_master_fixed <- merge(egenes_master, m2, by=c('gene_id','tissue'), all.x=T)
head(egenes_master_fixed)
egenes_master_fixed[is.na(D), D := 1.0]

save(egenes_master_fixed, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master_fixed.RData')