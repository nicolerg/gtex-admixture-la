library(data.table)
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

keep <- "/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt"
keep <- fread(keep, sep='\t', header=FALSE)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

args <- commandArgs(trailingOnly=TRUE)
cutoff <- as.numeric(args[1])

egenes_master <- egenes_master[pval_nominal < cutoff]
print(cutoff)

dflist <- list()
i <- 1 
for (t in tissues){

	cov <- paste0('/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/',t,'.v8.covariates.txt')

	egenes <- egenes_master[tissue==t]

	# get N admixed individuals 

	cov <- fread(cov, sep='\t', header=TRUE)
	subjs <- colnames(cov)[2:length(colnames(cov))]
	filt <- subjs[subjs %in% keep[,V1]]
	n_admixed <- length(filt)

	# get N global eGenes (1e-6)
	global <- egenes[method=='global']
	n_global <- length(unique(global[,gene_id]))

	# get N LAVA eGenes (1e-6)
	lava <- egenes[method=='LAVA']
	n_lava <- length(unique(lava[,gene_id]))

	# get N overlapping eGenes 
	n_overlapping <- length(unique(lava[gene_id %in% global[,gene_id], gene_id]))
	lava <- lava[gene_id %in% global[,gene_id]]
	global <- global[gene_id %in% lava[,gene_id]]

	# get N overlapping eGenes with different lead SNPs
	lava[,pair:=paste(gene_id, variant_id, sep=':')]
	global[,pair:=paste(gene_id, variant_id, sep=':')]

	repeats <- lava[pair %in% global[,pair], gene_id]

	n_diff_lead <- length(unique(lava[!(gene_id %in% repeats),gene_id]))

	print(t)
	print(paste('N admixed:', n_admixed))
	print(paste('N local eGenes:', n_lava))
	print(paste('N global eGenes:', n_global))
	print(paste('N overlapping:', n_overlapping))
	print(paste('N diff lead:',n_diff_lead))
	system("echo")

}

