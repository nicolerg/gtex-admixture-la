# compare GTEx eQTLs and LocalAA eQTLs

library(data.table)

master_data_dir = "/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/"
load(sprintf('%s/egenes_master-20200326.RData',master_data_dir))

head(egenes_master)

# egenes_master[LD == 1, overlapping_lead_variants := 0]
# write.table(egenes_master, file='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/egenes_master_supplement.tsv', sep='\t', col.names=T, row.names=F, quote=F)

egenes_master = egenes_master[pval_nominal_local < 1e-6]
head(egenes_master)
print(nrow(egenes_master))

gtex_in = '/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL'

i = 1
gtex_list = list()
for (t in unique(egenes_master[,tissue])){
  gtex = fread(cmd=sprintf('zcat %s/%s.v8.egenes.txt.gz', gtex_in, t))
  my_eqtl = unique(egenes_master[tissue==t,gene_id])
  gtex_sub = gtex[gene_id %in% my_eqtl]
  if(nrow(gtex_sub) != length(my_eqtl)){
    print(t)
  }
  gtex_sub[,tissue := t]
  gtex_list[[i]] = gtex_sub
  i = i + 1
}

gtex_local = rbindlist(gtex_list)
max(gtex_local[,pval_nominal])


master_data_dir = "/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/"
load(sprintf('%s/egenes_master-20200326.RData',master_data_dir))
egenes_master = egenes_master[pval_nominal_global < 1e-6]
head(egenes_master)
print(nrow(egenes_master))

gtex_in = '/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL'

i = 1
gtex_list = list()
for (t in unique(egenes_master[,tissue])){
  gtex = fread(cmd=sprintf('zcat %s/%s.v8.egenes.txt.gz', gtex_in, t))
  my_eqtl = unique(egenes_master[tissue==t,gene_id])
  gtex_sub = gtex[gene_id %in% my_eqtl]
  if(nrow(gtex_sub) != length(my_eqtl)){
    print(t)
  }
  gtex_sub[,tissue := t]
  gtex_list[[i]] = gtex_sub
  i = i + 1
}

gtex_global = rbindlist(gtex_list)
max(gtex_global[,pval_nominal])


gtex_both = merge(gtex_global, gtex_local, by=c('gene_id','tissue'), suffixes = c('.global','.local'), all = T)
nrow(gtex_both)

nrow(gtex_both[!complete.cases(gtex_both)])
table(gtex_both[,pval_nominal.global] < 1e-6, gtex_both[,pval_nominal.local] < 1e-6)

nrow(gtex_both[pval_nominal.local > 1e-6 | pval_nominal.global > 1e-6])

not_called_local = gtex_both[pval_nominal.local > 1e-6, .(tissue, gene_id)] # 270
not_called_global = gtex_both[pval_nominal.global > 1e-6, .(tissue, gene_id)] # 261 6938

not_called_global[,pair := paste(tissue, gene_id, ':')]
not_called_local[,pair := paste(tissue, gene_id, ':')]

table(not_called_local[,pair] %in% not_called_global[,pair])

nrow(gtex_both[pval_nominal.local > 1e-6])
nrow(gtex_both[pval_nominal.global > 1e-6])

