library(data.table)

# parse allpairs for custom colocalization with COLOC

# args <- commandArgs(trailingOnly=TRUE)
# tissue <- args[1]

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

# # get subset of genes we're interested in 

# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/diff_lead-1e-06.RData') # same eGene, different lead SNP 

# for (tissue in tissues){
# 	gene_list <- data.table(gene_id=unique(diff_lead[tissue==tissue & (p3 > 0.8 | p4 > 0.5),gene_id]))
# 	trait_counts <- data.table(table(diff_lead[tissue==tissue & (p3 > 0.8 | p4 > 0.5),trait]))
# 	write.table(trait_counts, paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/genelist/',tissue,'.traits.txt'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

# 	# filter allpairs 

# 	write.table(gene_list, paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/genelist/',tissue,'.genelist.txt'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
# 	system(sprintf('python ~/gtex-admix/scripts/LAVA/annotation/parse_allpairs_for_coloc_smr.py %s &',tissue))
# }

# quit()

# combine all SNPs into one file 

# system('cat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/af/*.snps.txt | sort | uniq > /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/af/allsnps.txt')
# snp_list <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/af/allsnps.txt' 

# get allele frequencies 

# geno_base <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/'

# for (chrom in 1:22){
# 	outprefix <- sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/af/chr%s.af',chrom)
# 	vcf <- sprintf('%sgtex.admixed.MAC10.phased.chr%s.recode.vcf.gz', geno_base, chrom)
# 	system(sprintf('vcftools --gzvcf %s --out %s --snps %s --recode &', vcf, outprefix, snp_list))
# 	system(sprintf('vcftools --gzvcf %s --out %s --snps %s --freq &', vcf, outprefix, snp_list))
# }

# add minor allele frequencies to allpairs files 

af_map <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/af/minor.af.tsv',sep='\t',header=TRUE)

for (tissue in tissues){

	local_out = sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/%s.filtered.local.allpairs.tsv.gz',tissue)
	global_out = sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/%s.filtered.global.allpairs.tsv.gz',tissue)

	allp_local <- fread(global_out, sep='\t', header=TRUE)
	allp_global <- fread(local_out, sep='\t', header=TRUE)

	allp_local <- merge(allp_local, af_map, by='variant_id')
	allp_global <- merge(allp_global, af_map, by='variant_id')

	write.table(allp_local, sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/%s.filtered.af.local.allpairs.tsv.gz',tissue), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
	write.table(allp_global, sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/smr-coloc/allpairs/%s.filtered.af.global.allpairs.tsv.gz',tissue), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
}

