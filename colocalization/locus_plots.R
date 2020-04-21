#!/bin/R

# make tables to generate LocusCompare-like plots 

library(data.table)
library(plyr)
library(ggplot2)

load('/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/master_coloc-for-plotting.RData')
# ID subset of loci 
m_sub = master_coloc[unique_hit == 1]
nrow(m_sub) # 31

table(m_sub[,ref_snp_global_finemap] == m_sub[,ref_snp_local_finemap])
table(m_sub[,ref_snp_global_coloc] == m_sub[,ref_snp_local_coloc])
table(m_sub[,ref_snp_global_coloc] == m_sub[,ref_snp_global_finemap])

m_sub[,global_better := ifelse(coloc_h4_global > coloc_h4_local, 1, 0)]

sites = m_sub[,.(gwas,gwas_trait,tissue,ref_snp_global_finemap,gene_name,gene_id,global_better)]
setnames(sites, 'ref_snp_global_finemap','ref_snp')

gwas_indir = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/gwas'
#imputed_UKB_20002_1463_self_reported_ulcerative_colitis.formatted.txt.gz
eqtl_indir = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/allpairs_difflead'
#Muscle_Skeletal.filtered.local.allpairs.sorted.tsv.gz
outdir='/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/locusplots'

plot_df <- list()
for (i in 1:nrow(sites)){
  
  print(i)
  
	pos <- as.integer(unlist(unname(strsplit(sites[i,ref_snp],'_')))[2])
	chr <- unlist(unname(strsplit(sites[i,ref_snp],'_')))[1]
	start <- pos - 1000000
	end <- pos + 1000000
	eqtl_locus <- paste0(chr,':',start,'-',end)
	gwas_locus <- paste0('chr',eqtl_locus)
	target_gene <- sites[i,gene_id]
	gene_label <- sites[i,gene_name]

	global_better <- sites[i,global_better]
	if(global_better == 1){
	  better_method = 'GlobalAA'
	}else{
	  better_method = 'LocalAA'
	}

	tissue <- sites[i,tissue]

	# get stats from Global eQTL file 
	geqtl <- sprintf('%s/%s.filtered.global.allpairs.sorted.tsv.gz',eqtl_indir,tissue)
	gout <- sprintf('%s/%s-%s-global.txt',outdir,gsub(':','_',eqtl_locus),tissue)
	system(sprintf('tabix %s %s > %s',geqtl,eqtl_locus,gout))

	# get stats form Local eQTL file
	lout <- sprintf('%s/%s-%s-local.txt',outdir,gsub(':','_',eqtl_locus),tissue)
	leqtl <- sprintf('%s/%s.filtered.local.allpairs.sorted.tsv.gz',eqtl_indir,tissue)
	system(sprintf('tabix %s %s > %s',leqtl,eqtl_locus,lout))

	# get stats from GWAS
	gwas_t <- sites[i,gwas]
	gwas_file <- sprintf('%s/%s.formatted.txt.gz',gwas_indir,gwas_t)
	gwasout <- sprintf('%s/%s-%s-global.txt',outdir,gsub(':','_',gwas_locus),gwas_t)
	system(sprintf('tabix %s %s > %s',gwas_file,gwas_locus,gwasout))

	gwas_ss <- fread(gwasout,sep='\t',header=FALSE)
	local_ss <- fread(lout,sep='\t',header=FALSE)
	global_ss <- fread(gout,sep='\t',header=FALSE)

	colnames(local_ss) <- c('chr','snp_pos','alt','ref','beta','se','pvalue','gene','effect_af')
	colnames(global_ss) <- c('chr','snp_pos','alt','ref','beta','se','pvalue','gene','effect_af')
	local_ss[,method := 'LocalAA']
	global_ss[,method := 'GlobalAA']
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

	#sub_diff <- diff_lead[tissue==tissue&gene==target_gene,.(method,variant_id,snp_pos,chr,alt,ref)]

	merged = m
	# merged <- merge(m, sub_diff, by=c('method','snp_pos','chr','ref','alt'),all.x=TRUE)
	# merged[,global_p4 := global_p4]
	# merged[,local_p4 := local_p4]

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

	merged[,description := sprintf('%s & %s: %s',gwas_t, tissue, gene_label)]
	merged[,description := sapply(description, fix_strings)]
	merged[,better_coloc := better_method]

	print(head(merged))

	plot_df[[i]] <- merged
	i <- i +1 
}

master <- data.table(rbindlist(plot_df))
save(master, file=sprintf('%s/master-locus_plots.RData',outdir))
