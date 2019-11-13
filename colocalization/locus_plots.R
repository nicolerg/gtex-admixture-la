#!/bin/R

# make tables to generate LocusCompare-like plots 

library(data.table)
library(plyr)
library(ggplot2)

# already filtered by pval and lead SNPs
load('coloc/master_coloc-20191030.RData') # colocalization results

gene_names <- fread('gene_name_map.tsv', sep='\t', header=TRUE) # ENSEMBL ID to gene name conversion

master_coloc[, gene_stable_id := sapply(gene_id, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)

master_coloc[,label_global:= ifelse(clpp_h4_local < 0.50 & clpp_h4_global > 0.5 & (clpp_h4_global - clpp_h4_local) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
master_coloc[,label_local:= ifelse(clpp_h4_local > 0.5 & clpp_h4_global < 0.5 & (clpp_h4_local - clpp_h4_global) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

global_better <- master_coloc[!is.na(label_global)]
local_better <- master_coloc[!is.na(label_local)]

# now get "seed SNP" of each eGene tested
load('merged/egenes_master-20191030.RData') # lead variants for all tested genes, called egenes_master

# filter by LD 
egenes = fread('plink-r2/20191030/evariant_pairs_test_ld.txt', sep='\t', header=F)
colnames(egenes) = c('gene_id','tissue','lead_local','lead_global')

global = egenes[,.(gene_id,tissue,lead_global)]
local = egenes[,.(gene_id,tissue,lead_local)]

colnames(global) = c('gene_id','tissue','variant_id')
global[,method := 'Global']
colnames(local) = c('gene_id','tissue','variant_id')
local[,method := 'Local']

diff_lead = rbindlist(list(global,local))

diff_lead[,snp_pos := as.integer(unname(unlist(strsplit(variant_id,'_')))[2])]
diff_lead[,chr := unname(unlist(strsplit(variant_id,'_')))[1]]
diff_lead[,alt := unname(unlist(strsplit(variant_id,'_')))[4]]
diff_lead[,ref := unname(unlist(strsplit(variant_id,'_')))[3]]
diff_lead[,gene := gene_id]
diff_lead[,gene_id := NULL]

sites <- data.table(rbind(global_better[,.(gwas_trait,tissue,ref_snp,gene_name,gene_id,clpp_h4_global,clpp_h4_local)],local_better[,.(gwas_trait,tissue,ref_snp,gene_name,gene_id,clpp_h4_global,clpp_h4_local)]))

plot_df <- list()
i <- 1
for (i in 1:nrow(sites)){
	pos <- as.integer(unlist(unname(strsplit(sites[i,ref_snp],'_')))[2])
	chr <- unlist(unname(strsplit(sites[i,ref_snp],'_')))[1]
	start <- pos - 500000
	end <- pos + 500000
	eqtl_locus <- paste0(chr,':',start,'-',end)
	gwas_locus <- paste0('chr',eqtl_locus)
	target_gene <- sites[i,gene_id]
	gene_label <- sites[i,gene_name]

	global_p4 <- sites[i,clpp_h4_global]
	local_p4 <- sites[i,clpp_h4_local]

	tissue <- sites[i,tissue]

	# get stats from Global eQTL file 
	geqtl <- sprintf('coloc/allpairs_difflead/%s.filtered.global.allpairs.sorted.tsv.gz',tissue)
	gout <- sprintf('coloc/locus_plots/%s-%s-global.txt',gsub(':','_',eqtl_locus),tissue)
	system(sprintf('tabix %s %s > %s',geqtl,eqtl_locus,gout))

	# get stats form Local eQTL file
	lout <- sprintf('coloc/locus_plots/%s-%s-local.txt',gsub(':','_',eqtl_locus),tissue)
	leqtl <- sprintf('coloc/allpairs_difflead/%s.filtered.local.allpairs.sorted.tsv.gz',tissue)
	system(sprintf('tabix %s %s > %s',leqtl,eqtl_locus,lout))

	# get stats from GWAS
	gwas <- sites[i,gwas_trait]
	gwas_file <- sprintf('coloc/gtex_gwas/coloc_input/coloc_imputed_%s.txt.gz',gwas)
	gwasout <- sprintf('coloc/locus_plots/%s-%s-global.txt',gsub(':','_',gwas_locus),gwas)
	system(sprintf('tabix %s %s > %s',gwas_file,gwas_locus,gwasout))

	gwas_ss <- fread(gwasout,sep='\t',header=FALSE)
	local_ss <- fread(lout,sep='\t',header=FALSE)
	global_ss <- fread(gout,sep='\t',header=FALSE)

	colnames(local_ss) <- c('chr','snp_pos','alt','ref','beta','se','pvalue','gene','effect_af')
	colnames(global_ss) <- c('chr','snp_pos','alt','ref','beta','se','pvalue','gene','effect_af')
	local_ss[,method := 'Local']
	global_ss[,method := 'Global']
	m_eqtl <- data.table(rbind(local_ss, global_ss))
	m_eqtl[,chr:=paste0('chr',chr)]
	m_eqtl <- m_eqtl[,.(chr, snp_pos, alt, ref, pvalue, gene, method)]
	m_eqtl <- m_eqtl[gene==target_gene]

	colnames(gwas_ss) <- c('chr','snp_pos','alt','ref','beta','se','pvalue','effect_af')
	gwas_ss <- gwas_ss[,.(chr, snp_pos, alt, ref, pvalue)]
	
	m <- merge(m_eqtl, gwas_ss, by=c('chr','snp_pos','alt','ref'),suffixes=c('_eqtl','_gwas'))

	m[,pvalue_gwas := -log10(pvalue_gwas)]
	m[,pvalue_eqtl := -log10(pvalue_eqtl)]
	m <- m[pvalue_eqtl > 1]

	sub_diff <- diff_lead[tissue==tissue&gene==target_gene,.(method,variant_id,snp_pos,chr,alt,ref)]

	merged <- merge(m, sub_diff, by=c('method','snp_pos','chr','ref','alt'),all.x=TRUE)
	merged[,global_p4 := global_p4]
	merged[,local_p4 := local_p4]

	# fix strings
	fix_strings = function(x){
		x = gsub('self_reported','SR',x)
		x = gsub('diagnosed_by_doctor','DBD',x)
		x = gsub('Skin_Not_Sun_Exposed_Suprapubic','NSE_Skin',x)
		x = gsub('Muscle_Skeletal','Skeletal_Muscle',x)
		x = gsub('Adipose_Subcutaenous','Subc._Adipose',x)
		x = gsub('_',' ',x)
		return(x)
	}

	merged[,description := sprintf('%s & %s: %s',gwas, tissue, gene_label)]
	merged[,description := sapply(description, fix_strings)]

	p4_lab <- as.character(sprintf('GlobalAA PP4 = %s\nLocalAA PP4 = %s',round(unique(global_p4),digits=2),round(unique(local_p4),digits=2)))

	merged[,p4_lab := p4_lab]

	print(head(merged))

	plot_df[[i]] <- merged
	i <- i +1 
}

master <- data.table(rbindlist(plot_df))
save(master, file='coloc/master-locus_plots.RData')
