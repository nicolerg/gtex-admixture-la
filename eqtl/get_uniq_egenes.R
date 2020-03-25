library(data.table)

# 25 March 2019
# updated 20 May 2019 
# get eQTLs that are unique to one method OR eGenes unique to one method OR eGenes with different lead SNPs 
# save(uniq_egene,file=paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/uniq_egene_',cutoff,'.RData'))
# save(uniq_eqtl,file=paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/uniq_eqtl_',cutoff,'.RData'))
# save(diff_lead,file=paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_',cutoff,'.RData'))

# first, make egenes_master.RData

args <- commandArgs(trailingOnly=TRUE)
base <- args[1] # /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed
outdir <- args[2] # /mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/merged

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic',
	'Whole_Blood')

# make egenes_master and allp_master

if(!file.exists(paste0(outdir,'/egenes_master.RData')) | !file.exists(paste0(outdir,'/allp_master.RData'))){

	allp_dt <- list()
	egenes_dt <- list()

	i <- 1

	for (t in tissues){

		basedir = paste0(base,'/',t,'/')
		filt_allpairs = paste0(basedir,t,'-LocalAA-GlobalAA-allpairs-filt-merged.tsv.gz')
		global_egenes = paste0(basedir,'gtex.admix.global.egenes.tied.txt.gz')
		lava_egenes = paste0(basedir,'gtex.admix.lava.egenes.tied.txt.gz')

		allp_filt = fread(cmd=sprintf("zcat %s",filt_allpairs), sep='\t', header=TRUE)
		global = fread(cmd=sprintf("zcat %s",global_egenes), sep='\t', header=TRUE)
		lava = fread(cmd=sprintf("zcat %s",lava_egenes), sep='\t', header=TRUE)

		# make a merged allpairs file
		allp_filt[,tissue:=t]

		# make a merged egenes file 
		global[,method:='global']
		global[,tissue:=t]
		lava[,method:='LAVA']
		lava[,tissue:=t]
		egenes = data.table(rbind(global, lava))

		allp_dt[[i]] <- allp_filt
		egenes_dt[[i]] <- egenes

		i <- i + 1

	}

	allp_master <- data.table(rbindlist(allp_dt))
	egenes_master <- data.table(rbindlist(egenes_dt))
	egenes_master[,pval_nominal := as.numeric(pval_nominal)]

	save(allp_master, file=paste0(outdir,'/allp_master.RData'))
	save(egenes_master, file=paste0(outdir,'/egenes_master.RData'))
} else {
	load(paste0(outdir,'/egenes_master.RData'))
}

# get_egenes <- function(cutoff){
# 	egenes_master <- egenes_master[pval_nominal < cutoff]

# 	global <- egenes_master[method=='global']
# 	lava <- egenes_master[method=='LAVA']

# 	global <- global[,.(gene_id, pval_nominal, tissue, method)]
# 	lava <- lava[,.(gene_id, pval_nominal, tissue, method)]

# 	uniq_df <- list()
# 	i <- 1
# 	for (t in tissues){
# 		g_uniq <- global[tissue==t & !(gene_id %in% lava[tissue==t, gene_id])]
# 		l_uniq <- lava[tissue==t & !(gene_id %in% global[tissue==t, gene_id])]
# 		g_uniq <- unique(g_uniq)
# 		l_uniq <- unique(l_uniq)
# 		m <- data.table(rbind(g_uniq,l_uniq))
# 		uniq_df[[i]] <- m
# 		i <- i+1
# 	}

# 	master_uniq_per_tissue <- data.table(rbindlist(uniq_df))
# 	master_uniq_genes <- data.frame(gene_id=unique(master_uniq_per_tissue[,gene_id]))

# 	#write.table(master_uniq_per_tissue, paste0(outdir,'/uniq_egenes_per_tissue_',cutoff,'.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
# 	write.table(master_uniq_genes, paste0(outdir,'/uniq_egenes_list_',cutoff,'.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
# }

#get_egenes(1e-4)
#get_egenes(1e-6)

filter_master <- function(cutoff){
	
	egenes_master <- egenes_master[pval_nominal < cutoff]

	uniq_egene <- list()
	uniq_eqtl <- list()
	diff_lead <- list()
	i <- 1
	for (t in tissues){

		egenes_master[,pair := paste0(gene_id, ':', variant_id)]

		global <- egenes_master[method=='global' & tissue == t]
		lava <- egenes_master[method=='LAVA' & tissue == t]

		# all unique eQTLs (including same eGene diff lead SNP)
		repeats <- lava[ pair %in% global[,pair] , gene_id ]
		lava <- lava[!(gene_id %in% repeats)]
		global <- global[!(gene_id %in% repeats)]
		#uniq_eqtl[[i]] <- data.table(rbind(global, lava)) # this doesn't test the case if one method's eSNP is a significant but not lead SNP in the other method 

		# now only same eGenes, different lead SNPs 
		l <- lava[gene_id %in% global[,gene_id]]
		g <- global[gene_id %in% lava[,gene_id]]
		diff_lead[[i]] <- data.table(rbind(g,l))

		# now only unique eGenes
		l <- lava[!(gene_id %in% global[,gene_id])]
		g <- global[!(gene_id %in% lava[,gene_id])]
		uniq_egene[[i]] <- data.table(rbind(g,l))

		i <- i + 1 

	}

	uniq_egene <- data.table(rbindlist(uniq_egene))
	#uniq_eqtl <- data.table(rbindlist(uniq_eqtl))
	diff_lead <- data.table(rbindlist(diff_lead))

	save(uniq_egene,file=paste0(outdir,'/uniq_egene_',cutoff,'.RData'))
	#save(uniq_eqtl,file=paste0(outdir,'/uniq_eqtl_',cutoff,'.RData'))
	save(diff_lead,file=paste0(outdir,'/diff_lead_',cutoff,'.RData'))

	# write out a file of same egenes, different lead SNP 
	same_gene_diff_lead <- data.table(gene_id=unique(diff_lead[,gene_id]))
	write.table(same_gene_diff_lead, file=paste0(outdir,'/same_egene_diff_lead_',cutoff,'.txt'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

}

filter_master(1e-4)
filter_master(1e-6)

