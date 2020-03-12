library(data.table)
library(ggplot2)

# see how much variance max Fst explains in coloc p4 

# add in continent-specific maxes

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/merged-all_no_max.RData')
# head(merged)
# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-04.RData') # eQTL calling

# diff_lead <- unique(diff_lead[,.(gene_id,tissue,method,variant_id)])
# diff_lead[,CHROM := sapply(variant_id, function(x) gsub('chr','',unname(unlist(strsplit(x,'_')))[1]))]
# diff_lead[,POS := sapply(variant_id, function(x) unname(unlist(strsplit(x,'_')))[2])]
# head(diff_lead)
# egene_fst <- merge(diff_lead, merged, by=c('CHROM','POS'))

# save(egene_fst,file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst-merged_no_max_diff_lead.RData')

# head(egene_fst)

# eur_cols <- colnames(egene_fst)[grepl('EUR',colnames(egene_fst))]
# afr_cols <- colnames(egene_fst)[grepl('AFR',colnames(egene_fst))]

# other_cols <- c('gene_id','tissue','method','variant_id')
# print('0')
# egene_fst[, (c(eur_cols,afr_cols)) := lapply(.SD, as.numeric), .SDcols = c(eur_cols,afr_cols)]
# print('0.1')
# egene_fst[,max_eur_fst := max(abs(.SD),na.rm=TRUE), by=c('gene_id','tissue','method'), .SDcols = eur_cols]
# print('0.2')
# egene_fst[,max_afr_fst := max(abs(.SD),na.rm=TRUE), by=c('gene_id','tissue','method'), .SDcols = afr_cols]
# print('1')
# save(egene_fst, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/egene_fst-max_per_superpop.RData')

# remove genes that don't have different lead SNPs
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master-20191030.RData')
egenes_master = egenes_master[,.(tissue, gene_id, LD, overlapping_lead_variants, pval_nominal_local, pval_nominal_global)]
egenes_master = egenes_master[pval_nominal_local < 1e-4 & pval_nominal_global < 1e-4]
egenes_master[,pval_diff := -log10(pval_nominal_local) + log10(pval_nominal_global)]
egenes_master = egenes_master[LD < 1 & overlapping_lead_variants == 0]
egenes_master[,LD := NULL]
egenes_master[,overlapping_lead_variants:=NULL]
# now it only includes "different lead variants"

## get FST info for each gene 

afr_eur_fst <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/high_local_fst/EUR_AFR.weir.fst',sep='\t',header=TRUE)
afr_eur_fst <- afr_eur_fst[WEIR_AND_COCKERHAM_FST == '-nan', WEIR_AND_COCKERHAM_FST := 0]
afr_eur_fst <- afr_eur_fst[WEIR_AND_COCKERHAM_FST != 'WEIR_AND_COCKERHAM_FST']
afr_eur_fst[,WEIR_AND_COCKERHAM_FST := as.numeric(WEIR_AND_COCKERHAM_FST)]

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/egene_fst-max_per_superpop.RData')
egene_fst <- egene_fst[,.(CHROM,POS,gene_id,tissue,method,max_eur_fst,max_afr_fst)]
egene_fst <- egene_fst[max_eur_fst == Inf | max_eur_fst == -Inf,  max_eur_fst := 0]
egene_fst <- egene_fst[max_afr_fst == Inf | max_afr_fst == -Inf,  max_eur_fst := 0]

# first, get max afr_eur_fst for each gene 
genemap = egene_fst[,.(CHROM,POS,gene_id,tissue,method)]
afr_eur_fst1 = merge(afr_eur_fst, genemap, by=c('CHROM','POS'))

print(nrow(genemap))
print(nrow(afr_eur_fst1))

afr_eur_fst2 = afr_eur_fst1[,list(WEIR_AND_COCKERHAM_FST = max(WEIR_AND_COCKERHAM_FST)), by=.(gene_id, tissue, method)]

# now, get max within-continent FST for each gene
egene_fst2 = egene_fst[,list(max_eur_fst=max(max_eur_fst),max_afr_fst=max(max_afr_fst)),by=.(gene_id,tissue,method)]

# now merge them, for a single set of values per gene 
merged <- merge(egene_fst2, afr_eur_fst2, by=c('gene_id','tissue','method'))

# make wider 
l <- merged[method=='LAVA']
g <- merged[method=='global']
l[,method:=NULL]
g[,method:=NULL]
merged <- merge(g, l, by=c('gene_id','tissue'),suffixes=c('_global','_local'))
# only 5499 left. is this is right set?

# try merging with egenes_master
m2 = merge(merged, egenes_master, by=c('gene_id','tissue'))
# 3694 out of 4019
m2[,max_within_continent_fst_local := mapply(max, max_eur_fst_local, max_afr_fst_local)]
m2[,max_within_continent_fst_global := mapply(max, max_eur_fst_global, max_afr_fst_global)]

# try merging with coloc. see how much I am missing 
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master_coloc-20191030.RData')
master_coloc[,LD := NULL]
#max_egene_fst <- egene_fst[,list(max_eur_fst_per_feature=max(max_eur_fst), max_afr_fst_per_feature=max(max_afr_fst)), by = .(method, gene_id, tissue)]

print(nrow(master_coloc))
master_coloc2 = merge(master_coloc, egenes_master[,.(gene_id,tissue)], by=c('gene_id','tissue'))
# 3920 out of 4019 tested 
# 57786 total tests 

m3 = merge(master_coloc2, m2, by=c('gene_id','tissue'))
# 53411 out of 57786 remaining

m3[,p4_diff := clpp_h4_local - clpp_h4_global]
m3[,within := max_within_continent_fst_local - max_within_continent_fst_global]
m3[,between := WEIR_AND_COCKERHAM_FST_local - WEIR_AND_COCKERHAM_FST_global]

# get max coloc per gene per tissue per method 
m4 = m3[,list(clpp_h4_global = max(clpp_h4_global),
	clpp_h4_local = max(clpp_h4_local),
	gwas_trait_local = which.max(clpp_h4_local),
	gwas_trait_global = which.max(clpp_h4_global)), by=.(gene_id, tissue, max_within_continent_fst_global, max_within_continent_fst_local,WEIR_AND_COCKERHAM_FST_global, WEIR_AND_COCKERHAM_FST_local, pval_nominal_global, pval_nominal_local)]

m4[,pval_nominal_global := -log10(pval_nominal_global)]
m4[,pval_nominal_local := -log10(pval_nominal_local)]

summary = summary(lm(data=m4, formula=pval_nominal_global ~ max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_global))
g = data.frame(summary$coefficients)
print(summary$r.squared)
g$variable <- rownames(g)
g = data.table(g)

summary = summary(lm(data=m4, formula=pval_nominal_local ~ max_within_continent_fst_local + WEIR_AND_COCKERHAM_FST_local))
l = data.frame(summary$coefficients)
print(summary$r.squared)
l$variable <- rownames(l)
l = data.table(l)

results = rbindlist(list(g,l))

colnames(results) = c('estimate','stderr','t','pvalue','variable')

print(results)

summary = summary(lm(data=m4, formula=clpp_h4_global ~ max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_global))
g = data.frame(summary$coefficients)
print(summary$r.squared)
g$variable <- rownames(g)
g = data.table(g)

summary = summary(lm(data=m4, formula=clpp_h4_local ~ max_within_continent_fst_local + WEIR_AND_COCKERHAM_FST_local))
l = data.frame(summary$coefficients)
print(summary$r.squared)
l$variable <- rownames(l)
l = data.table(l)

results = rbindlist(list(g,l))

colnames(results) = c('estimate','stderr','t','pvalue','variable')

print(results)

summary = summary(lm(data=m4, formula=clpp_h4_global ~ pval_nominal_global + max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_global))
g = data.frame(summary$coefficients)
print(summary$r.squared)
g$variable <- rownames(g)
g = data.table(g)

summary = summary(lm(data=m4, formula=clpp_h4_local ~ pval_nominal_local + max_within_continent_fst_local + WEIR_AND_COCKERHAM_FST_local))
l = data.frame(summary$coefficients)
print(summary$r.squared)
l$variable <- rownames(l)
l = data.table(l)

results = rbindlist(list(g,l))

colnames(results) = c('estimate','stderr','t','pvalue','variable')

print(results)

q()

write.table(results,'/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/merged-regression-results6.txt',sep='\t',col.names=T, row.names=F, quote=F)

# dtlistglobal <- list()
# dtlistlocal <- list()
# i <- 1

# for(filter in c(0,0.2,0.4,0.6,0.8)){

# 	print(filter)

# 	# global 
# 	g = m4[clpp_h4_global > filter]
# 	summary = summary(lm(data=g, formula=clpp_h4_global ~ max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_global))
# 	summ = data.frame(summary$coefficients)
# 	print(summary$r.squared)
# 	summ$variable <- rownames(summ)
# 	summ = data.table(summ)
# 	colnames(summ) = c('estimate','stderr','t','pvalue','variable')
# 	summ[,filter := filter]
# 	dtlistglobal[[i]] <- summ
	
# 	# local
# 	l = m4[clpp_h4_local > filter]
# 	summary = summary(lm(data=l, formula=clpp_h4_local ~ max_within_continent_fst_local + WEIR_AND_COCKERHAM_FST_local))
# 	summ = data.frame(summary$coefficients)
# 	print(summary$r.squared)
# 	summ$variable <- rownames(summ)
# 	summ = data.table(summ)
# 	colnames(summ) = c('estimate','stderr','t','pvalue','variable')
# 	summ[,filter := filter]
# 	dtlistlocal[[i]] <- summ

# 	i <- i + 1
# }

# global <- rbindlist(dtlistglobal)
# local <- rbindlist(dtlistlocal)

# global[,type := 'GlobalAA']
# local[,type := 'LocalAA']

# results = rbindlist(list(global, local))

# results

# write.table(results,'/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/merged-regression-results5.txt',sep='\t',col.names=T, row.names=F, quote=F)


# for(filter in c(0,0.1,0.2,0.3,0.4,0.5)){

# 	print(filter)

# 	m2 = m3[clpp_h4_local > filter | clpp_h4_global > filter]
	
# 	summ = data.frame(summary(lm(data=m2, formula=p4_diff ~ max_within_continent_fst_local + max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_local + WEIR_AND_COCKERHAM_FST_global))$coefficients)
# 	print(summary(lm(data=m2, formula=p4_diff ~ max_within_continent_fst_local + max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_local + WEIR_AND_COCKERHAM_FST_global))$r.squared)
# 	summ$variable <- rownames(summ)
# 	summ = data.table(summ)
# 	colnames(summ) = c('estimate','stderr','t','pvalue','variable')
# 	summ[,filter := filter]
# 	dtlist1[[i]] <- summ

# 	summ = data.frame(summary(lm(data=m2, formula=p4_diff ~ within + between))$coefficients)
# 	summ$variable <- rownames(summ)
# 	summ = data.table(summ)
# 	colnames(summ) = c('estimate','stderr','t','pvalue','variable')
# 	summ[,filter := filter]
# 	dtlist2[[i]] <- summ

# 	m4 = m3[clpp_h4_local > filter]
# 	summ = data.frame(summary(lm(data=m4, formula=clpp_h4_local ~ max_within_continent_fst_local + WEIR_AND_COCKERHAM_FST_local ))$coefficients)
# 	summ$variable <- rownames(summ)
# 	summ = data.table(summ)
# 	colnames(summ) = c('estimate','stderr','t','pvalue','variable')
# 	summ[,filter := filter]
# 	dtlistlocal[[i]] <- summ 	

# 	m5 = m3[clpp_h4_global > filter]
# 	summ = data.frame(summary(lm(data=m5, formula=clpp_h4_global ~ max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_global ))$coefficients)
# 	summ$variable <- rownames(summ)
# 	summ = data.table(summ)
# 	colnames(summ) = c('estimate','stderr','t','pvalue','variable')
# 	summ[,filter := filter]
# 	dtlistglobal[[i]] <- summ

# 	i <- i + 1
# }

# dt1 <- rbindlist(dtlist1)
# dt2 <- rbindlist(dtlist2)
# global <- rbindlist(dtlistglobal)
# local <- rbindlist(dtlistlocal)

# dt1[,type := 'original']
# dt2[,type := 'delta']
# global[,type := 'GlobalAA']
# local[,type := 'LocalAA']

# results = rbindlist(list(dt1, dt2, global, local))

# write.table(results,'/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/merged-regression-results3.txt',sep='\t',col.names=T, row.names=F, quote=F)

# # check what happens when I test pval_diff instead

# m4 = unique(m3[,.(gene_id,tissue,max_within_continent_fst_local, max_within_continent_fst_global, WEIR_AND_COCKERHAM_FST_local , WEIR_AND_COCKERHAM_FST_global, pval_diff)])
# print(nrow(m4))
# print(nrow(unique(m4[,.(gene_id,tissue)])))

# summary(lm(data=m4, formula=pval_diff ~ max_within_continent_fst_local + max_within_continent_fst_global + WEIR_AND_COCKERHAM_FST_local + WEIR_AND_COCKERHAM_FST_global))

