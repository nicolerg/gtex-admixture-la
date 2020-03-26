library(data.table)

master_data_dir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/merged'

# make a master file of eQTLs, AFR-EUR Fst, within-AFR Fst, within-EUR Fst, r2 with local ancestry
# match results by gene and tissue for each method 
# what about CADD, TF binding, coloc, het_fraction?

if(!file.exists(sprintf("%s/f_all-merged-collapsed-wide.RData",master_data_dir))){

	if(!file.exists(sprintf('%s/f_all-merged-tmp.RData',master_data_dir))){

		# eQTL results (top for each gene; no cutoff)
		load(sprintf("%s/egenes_master.RData",master_data_dir))

		# within-continent Fst - all SNPs
		writeLines('merge eGenes with within-continent Fst...')
		load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/merged-all_no_max.RData')
		egenes_master[,CHROM := sapply(variant_id, function(x) unname(unlist(strsplit(x, '_')))[1])]
		egenes_master[,POS := sapply(variant_id, function(x) as.numeric(unname(unlist(strsplit(x, '_')))[2]))]
		merged[,CHROM := paste0('chr',CHROM)]
		merged[,POS := as.numeric(POS)]
		f <- merge(egenes_master, merged, by=c('CHROM','POS'), all.x=T)

		print(head(f))

		# AFR-v-EUR Fst - all SNPs
		writeLines('merge result with EUR-AFR Fst...')
		afr_eur <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/high_local_fst/EUR_AFR.weir.fst',sep='\t',header=T)
		colnames(afr_eur) <- c('CHROM','POS','EUR_AFR_FST')
		afr_eur[,CHROM := paste0('chr',CHROM)]
		afr_eur[,POS := as.numeric(POS)]
		f <- merge(f, afr_eur, by=c('CHROM','POS'),all.x=T)

		print(head(f))
		
		# r2 with local ancestry 
		writeLines('merge result with R-sq with LA...')
		r2 <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/REVISIONS/REVISED_all_lead_snp_r2.txt.gz',sep='\t',header=T)
		colnames(r2) <- c('variant_id','rsq_genotype_local_ancestry')
		f <- merge(f, r2, by='variant_id',all.x=T)

		print(head(f))

		# allele frequencies (as determined by 1000 Genomes populations)
		writeLines('merge result with 1000 Genomes superpopulation AFs...')
		af <- fread(cmd='zcat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/af/gtex.admixed.MAC10.snps.af.txt.gz',sep='\t',header=T)
		af <- af[,.(EAS_AF,EUR_AF,AFR_AF,variant_id)]
		f <- merge(f, af, by=c('variant_id'))
		
		print(head(f))

		# CADD
		writeLines('merge result with CADD...')
		cadd <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cadd/gtex.admixed.MAC10.all.snps_cadd.tsv',sep='\t',header=F)
		colnames(cadd) <- c('variant_id','CADD_score')
		f <- merge(f, cadd, by=c('variant_id'), all.x=T)
		
		print(head(f))

		# TF binding
		writeLines('merge result with ENCODE TF binding...')
		# gastroc
		gastroc <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/gastrocnemius_medialis/gastrocnemius_medialis.tf.binding.tsv',sep='\t',header=F)
		head(gastroc)
		gastroc[,tissue := 'Muscle_Skeletal']
		gastroc <- gastroc[,.(V3, V7, tissue)]
		# adipose
		adipose <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/subcutaneous_adipose_tissue/subcutaneous_adipose_tissue.tf.binding.tsv',sep='\t',header=F)
		adipose[,tissue := 'Adipose_Subcutaneous']
		adipose <- adipose[,.(V3, V7, tissue)]
		# nerve
		nerve <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/tibial_nerve/tibial_nerve.tf.binding.tsv',sep='\t',header=F)
		nerve[,tissue := 'Nerve_Tibial']
		nerve <- nerve[,.(V3, V7, tissue)]
		# artery 
		artery <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/tibial_artery/tibial_artery.tf.binding.tsv',sep='\t',header=F)
		artery[,tissue := 'Artery_Tibial']
		artery <- artery[,.(V3, V7, tissue)]
		# lung 
		lung <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/lung/lung.tf.binding.tsv',sep='\t',header=F)
		lung[,tissue := 'Lung']
		lung <- lung[,.(V3, V7, tissue)]
		tf_binding <- data.table(rbindlist(list(gastroc, adipose, nerve, artery, lung)))
		colnames(tf_binding) <- c('variant_id','tf_binding_annotation','tissue')
		tf_binding <- unique(tf_binding)

		# no TF binding for blood
		f <- merge(f, tf_binding, by=c('variant_id','tissue'), all.x=T)
		
		print(head(f))

		# # het fraction 
		# het <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/het_fraction/diff-lead-pairs-by-tissue-1e04-het-frac.tsv',sep='\t',header=T)
		# het <- het[,.(variant_id, tissue, het_frac)]
		# colnames(het) <- c('variant_id','tissue','frac_heterozygous_genotypes')
		# f <- merge(f, het, by=c('variant_id','tissue'),all.x=T)
		
		# print(head(f))

		# convert numeric columns to numeric

		num_cols <- c('AFR_GWD-AFR_ESN',
			'AFR_LWK-AFR_ESN',
			'AFR_LWK-AFR_GWD',
			'AFR_MSL-AFR_ESN',
			'AFR_MSL-AFR_GWD',
			'AFR_MSL-AFR_LWK',
			'AFR_YRI-AFR_ESN',
			'AFR_YRI-AFR_GWD',
			'AFR_YRI-AFR_LWK',
			'AFR_YRI-AFR_MSL',
			'EUR_FIN-EUR_CEU',
			'EUR_GBR-EUR_CEU',
			'EUR_GBR-EUR_FIN',
			'EUR_IBS-EUR_CEU',
			'EUR_IBS-EUR_FIN',
			'EUR_IBS-EUR_GBR',
			'EUR_TSI-EUR_CEU',
			'EUR_TSI-EUR_FIN',
			'EUR_TSI-EUR_GBR',
			'EUR_TSI-EUR_IBS',
			'EUR_AFR_FST',
			'rsq_genotype_local_ancestry')
			#'frac_heterozygous_genotypes')

		f[,(num_cols) := lapply(.SD, as.numeric), .SDcols=num_cols]
		print(head(f))
		f <- unique(f)
		save(f, file=sprintf("%s/f_all-merged-tmp.RData",master_data_dir))

	} else {
		load(sprintf("%s/f_all-merged-tmp.RData",master_data_dir))
	}

	if(!file.exists(sprintf("%s/f_all-merged-collapsed.RData",master_data_dir))){

		# collapse values down to one line per gene/tissue pair

		f[,c('CHROM','POS'):=NULL]
		f_collapsed <- f[,list(lead_variants=paste(variant_id, collapse=';'),
								`AFR_GWD-AFR_ESN`=mean(`AFR_GWD-AFR_ESN`,na.rm=T),
								`AFR_LWK-AFR_ESN`=mean(`AFR_LWK-AFR_ESN`,na.rm=T),
								`AFR_LWK-AFR_GWD`=mean(`AFR_LWK-AFR_GWD`,na.rm=T),
								`AFR_MSL-AFR_ESN`=mean(`AFR_MSL-AFR_ESN`,na.rm=T),
								`AFR_MSL-AFR_GWD`=mean(`AFR_MSL-AFR_GWD`,na.rm=T),
								`AFR_MSL-AFR_LWK`=mean(`AFR_MSL-AFR_LWK`,na.rm=T),
								`AFR_YRI-AFR_ESN`=mean(`AFR_YRI-AFR_ESN`,na.rm=T),
								`AFR_YRI-AFR_GWD`=mean(`AFR_YRI-AFR_GWD`,na.rm=T),
								`AFR_YRI-AFR_LWK`=mean(`AFR_YRI-AFR_LWK`,na.rm=T),
								`AFR_YRI-AFR_MSL`=mean(`AFR_YRI-AFR_MSL`,na.rm=T),
								`EUR_FIN-EUR_CEU`=mean(`EUR_FIN-EUR_CEU`,na.rm=T),
								`EUR_GBR-EUR_CEU`=mean(`EUR_GBR-EUR_CEU`,na.rm=T),
								`EUR_GBR-EUR_FIN`=mean(`EUR_GBR-EUR_FIN`,na.rm=T),
								`EUR_IBS-EUR_CEU`=mean(`EUR_IBS-EUR_CEU`,na.rm=T),
								`EUR_IBS-EUR_FIN`=mean(`EUR_IBS-EUR_FIN`,na.rm=T),
								`EUR_IBS-EUR_GBR`=mean(`EUR_IBS-EUR_GBR`,na.rm=T),
								`EUR_TSI-EUR_CEU`=mean(`EUR_TSI-EUR_CEU`,na.rm=T),
								`EUR_TSI-EUR_FIN`=mean(`EUR_TSI-EUR_FIN`,na.rm=T),
								`EUR_TSI-EUR_GBR`=mean(`EUR_TSI-EUR_GBR`,na.rm=T),
								`EUR_TSI-EUR_IBS`=mean(`EUR_TSI-EUR_IBS`,na.rm=T),
								EUR_AFR_FST=mean(EUR_AFR_FST,na.rm=T),
								rsq_genotype_local_ancestry=mean(rsq_genotype_local_ancestry,na.rm=T),
								TG_EAS_AF=mean(EAS_AF, na.rm=T),
								TG_EUR_AF=mean(EUR_AF, na.rm=T),
								TG_AFR_AF=mean(AFR_AF, na.rm=T),
								max_CADD_score=max(CADD_score, na.rm=T),
								which_max_CADD_score=variant_id[which.max(CADD_score)],
								min_CADD_score=min(CADD_score, na.rm=T),
								which_min_CADD_score=variant_id[which.min(CADD_score)]),
								#tf_binding_annotation=ifelse(1 %in% tf_binding_annotation, 1, 0),
								#frac_heterozygous_genotypes=mean(frac_heterozygous_genotypes, na.rm=T)),
			by=c('tissue','gene_id','count','pval_nominal','slope','method')]
		print(head(f_collapsed))
		f_collapsed <- unique(f_collapsed)

		# get max within-continent Fst
		afr_cols <- c('AFR_GWD-AFR_ESN',
			'AFR_LWK-AFR_ESN',
			'AFR_LWK-AFR_GWD',
			'AFR_MSL-AFR_ESN',
			'AFR_MSL-AFR_GWD',
			'AFR_MSL-AFR_LWK',
			'AFR_YRI-AFR_ESN',
			'AFR_YRI-AFR_GWD',
			'AFR_YRI-AFR_LWK',
			'AFR_YRI-AFR_MSL')

		eur_cols <- c('EUR_FIN-EUR_CEU',
			'EUR_GBR-EUR_CEU',
			'EUR_GBR-EUR_FIN',
			'EUR_IBS-EUR_CEU',
			'EUR_IBS-EUR_FIN',
			'EUR_IBS-EUR_GBR',
			'EUR_TSI-EUR_CEU',
			'EUR_TSI-EUR_FIN',
			'EUR_TSI-EUR_GBR',
			'EUR_TSI-EUR_IBS')

		f_eur <- f_collapsed[,list(max_EUR_fst = max(abs(.SD),na.rm=TRUE),
			which_max_EUR_fst = colnames(.SD)[max.col(abs(.SD))]), by=c('gene_id','tissue','method'), .SDcols = eur_cols]

		f_afr <- f_collapsed[,list(max_AFR_fst = max(abs(.SD),na.rm=TRUE),
			which_max_AFR_fst = colnames(.SD)[max.col(abs(.SD))]), by=c('gene_id','tissue','method'), .SDcols = afr_cols]

		fst <- merge(f_eur, f_afr, by=c('gene_id','tissue','method'))
		f_collapsed <- merge(f_collapsed, fst, by=c('gene_id','tissue','method'))

		save(f_collapsed, file=sprintf("%s/f_all-merged-collapsed.RData",master_data_dir))
	} else {
		load(sprintf("%s/f_all-merged-collapsed.RData",master_data_dir))
	}

	# get one line per gene
	g <- f_collapsed[method=='global']
	g[,method := NULL]
	l <- f_collapsed[method=='LAVA']
	l[,method := NULL]

	f <- merge(g, l, by=c('tissue','gene_id'), suffixes=c('.global','.local'))
	head(f)
	save(f, file=sprintf("%s/f_all-merged-collapsed-wide.RData",master_data_dir))
} else {
	load(sprintf("%s/f_all-merged-collapsed-wide.RData",master_data_dir))
}

nrow(f)

# merge with D of SNP pairs, when possible 
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/merged/egenes_master-20200326.RData')
egenes_master = egenes_master[,.(gene_id,tissue,LD,overlapping_lead_variants)]
colnames(egenes_master) = c('gene_id','tissue','lead_variants_LD','overlapping_lead_variants')

la_summary <- merge(f, egenes_master, by=c('gene_id','tissue'))
print(nrow(f))
print(nrow(la_summary))

save(la_summary, file="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/la_summary-all-annotations-merged.RData")

writeLines(colnames(la_summary))
# gene_id
# tissue
# count.global
# pval_nominal.global
# slope.global
# lead_variants.global
# AFR_GWD-AFR_ESN.global
# AFR_LWK-AFR_ESN.global
# AFR_LWK-AFR_GWD.global
# AFR_MSL-AFR_ESN.global
# AFR_MSL-AFR_GWD.global
# AFR_MSL-AFR_LWK.global
# AFR_YRI-AFR_ESN.global
# AFR_YRI-AFR_GWD.global
# AFR_YRI-AFR_LWK.global
# AFR_YRI-AFR_MSL.global
# EUR_FIN-EUR_CEU.global
# EUR_GBR-EUR_CEU.global
# EUR_GBR-EUR_FIN.global
# EUR_IBS-EUR_CEU.global
# EUR_IBS-EUR_FIN.global
# EUR_IBS-EUR_GBR.global
# EUR_TSI-EUR_CEU.global
# EUR_TSI-EUR_FIN.global
# EUR_TSI-EUR_GBR.global
# EUR_TSI-EUR_IBS.global
# EUR_AFR_FST.global
# rsq_genotype_local_ancestry.global
# TG_EAS_AF.global
# TG_EUR_AF.global
# TG_AFR_AF.global
# max_CADD_score.global
# which_max_CADD_score.global
# min_CADD_score.global
# which_min_CADD_score.global
# tf_binding_annotation.global
# max_EUR_fst.global
# which_max_EUR_fst.global
# max_AFR_fst.global
# which_max_AFR_fst.global
# count.local
# pval_nominal.local
# slope.local
# lead_variants.local
# AFR_GWD-AFR_ESN.local
# AFR_LWK-AFR_ESN.local
# AFR_LWK-AFR_GWD.local
# AFR_MSL-AFR_ESN.local
# AFR_MSL-AFR_GWD.local
# AFR_MSL-AFR_LWK.local
# AFR_YRI-AFR_ESN.local
# AFR_YRI-AFR_GWD.local
# AFR_YRI-AFR_LWK.local
# AFR_YRI-AFR_MSL.local
# EUR_FIN-EUR_CEU.local
# EUR_GBR-EUR_CEU.local
# EUR_GBR-EUR_FIN.local
# EUR_IBS-EUR_CEU.local
# EUR_IBS-EUR_FIN.local
# EUR_IBS-EUR_GBR.local
# EUR_TSI-EUR_CEU.local
# EUR_TSI-EUR_FIN.local
# EUR_TSI-EUR_GBR.local
# EUR_TSI-EUR_IBS.local
# EUR_AFR_FST.local
# rsq_genotype_local_ancestry.local
# TG_EAS_AF.local
# TG_EUR_AF.local
# TG_AFR_AF.local
# max_CADD_score.local
# which_max_CADD_score.local
# min_CADD_score.local
# which_min_CADD_score.local
# tf_binding_annotation.local
# max_EUR_fst.local
# which_max_EUR_fst.local
# max_AFR_fst.local
# which_max_AFR_fst.local
# lead_variants_LD
# overlapping_lead_variants