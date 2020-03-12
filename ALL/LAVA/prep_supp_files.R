library(data.table)

# prepare supplementary files

tissues <- c('Nerve_Tibial',
	'Artery_Tibial',
	'Adipose_Subcutaneous',
	'Muscle_Skeletal',
	'Skin_Not_Sun_Exposed_Suprapubic',
	'Lung'
	)

# egenes_master
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master-20191030.RData')
write.table(egenes_master, '~/gtex-admix/supplementary/LocalAA_GlobalAA_all_lead_variants.txt', sep='\t', col.names=T, row.names=F, quote=F)

# variance explained by LA and GA 
i <-1 
dt_list <- list()
for (t in tissues){
	dt <- fread(sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/%s_var_explained_by_ancestry.tsv',t),sep='\t',header=T)
	dt[,tissue := t]
	dt_list[[i]] <- dt
	i <- i + 1
}

merged <- rbindlist(dt_list)
write.table(merged, '~/gtex-admix/supplementary/LA_GA_VE_gene_expression.txt', sep='\t', col.names=T, row.names=F, quote=F)

# master_coloc
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master_coloc-20191030.RData')
master_coloc[,LD := NULL]
master_coloc[,overlapping_lead_variants := NULL]
write.table(master_coloc, '~/gtex-admix/supplementary/LocalAA_GlobalAA_colocalizations.txt', sep='\t', col.names=T, row.names=F, quote=F)
