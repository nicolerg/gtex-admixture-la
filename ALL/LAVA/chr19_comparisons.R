library(data.table)
library(ggplot2)

# 29 April 2019
# updated 20 May 2019

# compare Local results to an alternative method. either:
# 1. local + global in only admixed 
#	- chr19
#	- Muscle_Skeletal
#	- Skin_Not_Sun_Exposed_Suprapubic
# 2. local in all subjs (chr 19)
# 	- chr19
#	- Skin_Not_Sun_Exposed_Suprapubic (only this tissue for now)
# 3. local GTEx v8 peers vs. local recalculated peers in admixed only 

# # first, let's read in the chr19 standard results 

# my_allpairs_skin <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Skin_Not_Sun_Exposed_Suprapubic/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-merged-allpairs-chr19.tsv.gz'
# my_allpairs_muscle <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Muscle_Skeletal/Muscle_Skeletal-LAVA-global-merged-allpairs-chr19.tsv.gz'

# gtex_allpairs_skin <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/gtex_v8/Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz'
# gtex_allpairs_muscle <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/gtex_v8/Muscle_Skeletal.allpairs.txt.gz'

# allanc_allpairs_muscle <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/both_anc/Muscle_Skeletal-LAVA-global-merged-allpairs-chr19.tsv.gz'
# allanc_allpairs_skin <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/both_anc/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-merged-allpairs-chr19.tsv.gz'

# allsubj_allpairs_skin <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/all_ind/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-merged-allpairs-chr19.tsv.gz'
# # no muscle for all subjects with local ancestry yet

# what do I want to compare?
# 1. correlation of p-values and effect sizes for allpairs
# 2. if eGenes change (1e-6)
# 3. if lead SNPs change (1e-6)


# get lead SNPs for allpairs files
# output of chr19_combined_egene.py

compare_results <- function(my_method, other_method, cutoff, dt){
	## my method and all anc
	sub <- dt[,.(variant_id, gene_id, get(paste0('pval_nominal_',my_method)), get(paste0('slope_',my_method)), get(paste0('pval_nominal_',other_method)), get(paste0('slope_',other_method)))]
	colnames(sub) <- c('variant_id','gene_id',paste0('pval_nominal_',my_method),paste0('slope_',my_method),paste0('pval_nominal_',other_method),paste0('slope_',other_method))

	print(sprintf('n eGenes tested: %s', length(unique(sub[,gene_id]))))
	print(cutoff)
	sub <- sub[!is.na(get(paste0('pval_nominal_',my_method))) | !is.na(get(paste0('pval_nominal_',other_method)))]

	# filter down to 1e-6
	sub <- sub[ ( get(paste0('pval_nominal_',my_method)) < cutoff | is.na(get(paste0('pval_nominal_',my_method))) ) & (get(paste0('pval_nominal_',other_method)) < cutoff | is.na(get(paste0('pval_nominal_',other_method))) ) ]
	sub <- sub[order(gene_id)]

	#1. if eGenes change b/w my method and allanc
	my_method_egenes = unique(sub[!(is.na(get(paste0('pval_nominal_',my_method)))), gene_id])
	allanc_egenes = unique(sub[!(is.na(get(paste0('pval_nominal_',other_method)))), gene_id])

	egene_mymethod_only = length(setdiff(my_method_egenes, allanc_egenes)) # in x but not in y 
	egene_allanc_only = length(setdiff(allanc_egenes, my_method_egenes))

	# same egenes
	same_egenes <- sub[gene_id %in% intersect(my_method_egenes, allanc_egenes)]
	same_egene_diff_lead = length(same_egenes[duplicated(gene_id),gene_id])# duplicate gene_ids
	same_egene_same_lead = nrow(same_egenes[complete.cases(same_egenes)])# complete rows

	# if (my_method == 'allsubj_global' & other_method == 'gtex'){
	# 	x <- same_egenes[duplicated(gene_id)]
	# 	x <- same_egenes[order(gene_id)]
	# 	print(x)
	# }
	
	print(sprintf('n eGenes: %s', length(unique(sub[,gene_id]))))
	print(sprintf('n same eGenes: %s', length(unique(same_egenes[,gene_id]))))
	print(sprintf('n same eGenes diff lead: %s', same_egene_diff_lead))
	print(sprintf('n same eGene same lead: %s', same_egene_same_lead))
	print(sprintf('n unique eGenes my method: %s', egene_mymethod_only))
	print(sprintf('n unique eGenes other method: %s', egene_allanc_only))
	system('echo')

	counts <- data.table(N_diff_lead_snp = same_egene_diff_lead,
		N_same_lead_snp = same_egene_same_lead,
		N_method1_only = egene_mymethod_only,
		N_method2_only = egene_allanc_only,
		total_method1 = length(my_method_egenes),
		total_method2 = length(allanc_egenes))

	return(counts)

}

format_for_plotting <- function(dt){
	m <- dt[,.(tissue, N_diff_lead_snp, N_same_lead_snp, N_method1_only, N_method2_only, pval_cutoff)]

	m[,method1_y := N_method2_only + N_diff_lead_snp + N_same_lead_snp + 0.5*N_method1_only]
	m[,method2_y := 0.5*N_method2_only]

	m <- melt(m, id=c('tissue','method1_y','method2_y','pval_cutoff'), variable.name='type', value.name='count')
	m[,type := factor(type, levels = c('N_method1_only', 'N_same_lead_snp', 'N_diff_lead_snp' ,'N_method2_only'))]

	m[type == 'N_method2_only', method2_label := count]
	m[type == 'N_method1_only', method1_label := count]

	return(m)
}

plot_comparison <- function(method1, method2, dt){
	tissues <- c('Muscle_Skeletal',
		'Skin_Not_Sun_Exposed_Suprapubic')

	shortlab <- c(Muscle_Skeletal='Skeletal\nmuscle',
		Skin_Not_Sun_Exposed_Suprapubic='Skin (NSE)')

	labels <- c("Skeletal muscle"='Skeletal\nmuscle',
		"Skin (not sun-exposed)"="Skin (not\nsun-exposed)")

	dt[,pval_cutoff := factor(pval_cutoff)]

	g <- ggplot(dt) +
		geom_bar(aes(x=tissue, y=count, fill=type), stat='identity', colour='black') +
		geom_text(aes(x=tissue, y=method2_y, label=method2_label)) +
		geom_text(aes(x=tissue, y=method1_y, label=method1_label)) +
		scale_fill_manual(values=c(N_diff_lead_snp='#9999ff',
		                          N_same_lead_snp='gray',
		                          N_method2_only='#8DD3C7',
		                          N_method1_only='#FDB462'),
		                 labels=c(N_diff_lead_snp='  Same eGene, different lead SNP',
		                          N_same_lead_snp='  Same eGene, same lead SNP',
		                          N_method2_only=sprintf('  eGene in %s only',method2),
		                          N_method1_only=sprintf('  eGene in %s only',method1)),
		                 breaks=c('N_method1_only', 'N_same_lead_snp', 'N_diff_lead_snp' ,'N_method2_only')) +
		theme_bw() +
		theme(axis.title.x=element_blank(),
		   axis.text.x=element_text(colour='black'),
		   legend.title=element_blank(),
		   panel.grid.minor=element_blank(),
		   panel.border = element_rect(colour = "black", fill=NA, size=1),
		   legend.background = element_rect(color = "black", size = 0.5, linetype = "solid")) +
		scale_x_discrete(labels=shortlab) +
		labs(y='N eGenes') +
		facet_wrap(~ pval_cutoff, ncol=2)

	return(g)
}

ancestry_comparisons <- function(){
	muscle_compare <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/chr19-merged-all-comparisons-tied-egenes-Muscle_Skeletal.tsv.gz'
	skin_compare <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/chr19-merged-all-comparisons-tied-egenes-Skin_Not_Sun_Exposed_Suprapubic.tsv.gz'
	outdir <- '~/gtex-admix/plots/final_figures/'

	allanc_reproduce_plot <- list()
	i<-1
	gtex_reproduce_plot <- list()
	allsub_reproduce_plot <- list()
	j<-1

	for (file in c(skin_compare, muscle_compare)){

		dt <- fread(cmd=sprintf('zcat %s',file), sep='\t', header=TRUE)

		print(basename(file))
		
		for (cutoff in c(1e-6, 1e-5)){

			# compare my method (local) to all allanc
			print(sprintf('%s vs. %s','admixed_local', 'admixed_allanc'))
			allanc <- compare_results('admixed_local','admixed_allanc',cutoff,dt)
			allanc[,tissue := ifelse(file==skin_compare, 'Skin_Not_Sun_Exposed_Suprapubic','Muscle_Skeletal')]
			allanc[,pval_cutoff := cutoff]
			allanc_reproduce_plot[[i]] <- allanc
			i <- i + 1

			if (file == skin_compare){

				# compare my method to all subj (skin only)
				print(sprintf('%s vs. %s','admixed_local', 'allsubj_local'))
				allsub <- compare_results('admixed_local','allsubj_local',cutoff,dt)
				allsub[,tissue := 'Skin_Not_Sun_Exposed_Suprapubic']
				allsub[,pval_cutoff := cutoff]
				allsub_reproduce_plot[[j]] <- allsub 

				# compare my global results to GTEx v8 global results
				print(sprintf('%s vs. %s','admixed_global', 'gtex'))
				gtex_reproduce <- compare_results('allsubj_global','gtex',cutoff,dt)
				gtex_reproduce[,tissue := 'Skin_Not_Sun_Exposed_Suprapubic']
				gtex_reproduce[,pval_cutoff := cutoff]
				gtex_reproduce_plot[[j]] <- gtex_reproduce

				j <- j + 1

			}
		}


		if (file == skin_compare){

			gtex_plot <- data.table(rbindlist(gtex_reproduce_plot))
			allsub_plot <- data.table(rbindlist(allsub_reproduce_plot)) 

			gtex_final <- format_for_plotting(gtex_plot)
			allsub_final <- format_for_plotting(allsub_plot)

			g2 <- plot_comparison('all_subjects_global','gtex_v8',gtex_final)
			g3 <- plot_comparison('admixed_local','all_subjects_local',allsub_final)

			pdf(sprintf('%s/Skin-admixed-local-v-all-subjects-local.pdf',outdir),width=8, height=4)
			print(g3)
			dev.off()

			pdf(sprintf('%s/Skin-all-subjects-global-v-gtex-v8.pdf',outdir),width=8, height=4)
			print(g2)
			dev.off()
		}

	}

	allanc_plot <- data.table(rbindlist(allanc_reproduce_plot))
	allanc_final <- format_for_plotting(allanc_plot)

	g1 <- plot_comparison('admixed_local','admixed_all_ancestry',allanc_final)

	pdf(sprintf('%s/admixed-local-v-admixed-all-ancestry.pdf',outdir),width=8, height=4)
	print(g1)
	dev.off()
}



peer_comparisons <- function(){

	file <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/chr19-merged-newpeers-tied-egenes-Skin_Not_Sun_Exposed_Suprapubic.tsv.gz'
	outdir <- '~/gtex-admix/plots/final_figures/'

	peer_reproduce_plot <- list()
	i<-1

	dt <- fread(cmd=sprintf('zcat %s',file), sep='\t', header=TRUE)
	print(basename(file))
	
	for (cutoff in c(1e-6, 1e-5)){

		# compare old peers to new peers in Local
		print(sprintf('%s vs. %s','admixed_local', 'newpeers_local'))
		peer <- compare_results('admixed_local','newpeers_local',cutoff,dt)
		peer[,method := 'Local']
		peer[,pval_cutoff := cutoff]
		peer_reproduce_plot[[i]] <- peer
		i <- i + 1

		# compare old peers to new peers in Global
		print(sprintf('%s vs. %s','admixed_global', 'newpeers_global'))
		peer <- compare_results('admixed_global','newpeers_global',cutoff,dt)
		peer[,method := 'Global']
		peer[,pval_cutoff := cutoff]
		peer_reproduce_plot[[i]] <- peer
		i <- i + 1

	}

	peer_plot <- data.table(rbindlist(peer_reproduce_plot))

	m <- peer_plot[,.(method, N_diff_lead_snp, N_same_lead_snp, N_method1_only, N_method2_only, pval_cutoff)]

	m[,method1_y := N_method2_only + N_diff_lead_snp + N_same_lead_snp + 0.5*N_method1_only]
	m[,method2_y := 0.5*N_method2_only]

	m <- melt(m, id=c('method','method1_y','method2_y','pval_cutoff'), variable.name='type', value.name='count')
	m[,type := factor(type, levels = c('N_method1_only', 'N_same_lead_snp', 'N_diff_lead_snp' ,'N_method2_only'))]

	m[type == 'N_method2_only', method2_label := count]
	m[type == 'N_method1_only', method1_label := count]

	print(head(m))

	method1 <- 'GTEx v8 PEER factors'
	method2 <- 'Admixed PEER factors'

	m[,pval_cutoff := factor(pval_cutoff)]

	g <- ggplot(m) +
		geom_bar(aes(x=method, y=count, fill=type), stat='identity', colour='black') +
		geom_text(aes(x=method, y=method2_y, label=method2_label)) +
		geom_text(aes(x=method, y=method1_y, label=method1_label)) +
		scale_fill_manual(values=c(N_diff_lead_snp='#9999ff',
		                          N_same_lead_snp='gray',
		                          N_method2_only='#8DD3C7',
		                          N_method1_only='#FDB462'),
		                 labels=c(N_diff_lead_snp='  Same eGene, different lead SNP',
		                          N_same_lead_snp='  Same eGene, same lead SNP',
		                          N_method2_only=sprintf('  eGene in %s only',method2),
		                          N_method1_only=sprintf('  eGene in %s only',method1)),
		                 breaks=c('N_method1_only', 'N_same_lead_snp', 'N_diff_lead_snp' ,'N_method2_only')) +
		theme_bw() +
		theme(axis.title.x=element_blank(),
		   axis.text.x=element_text(colour='black'),
		   legend.title=element_blank(),
		   panel.grid.minor=element_blank(),
		   panel.border = element_rect(colour = "black", fill=NA, size=1),
		   legend.background = element_rect(color = "black", size = 0.5, linetype = "solid")) +
		labs(y='N eGenes') +
		facet_wrap(~ pval_cutoff, ncol=2)

	pdf(sprintf('%s/oldpeer-vs-newpeer-admixed-only.pdf',outdir),width=8, height=4)
	print(g)
	dev.off()
}
peer_comparisons()

