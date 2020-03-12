library(data.table)

# run HOMER on regions surrounding method-specific variants

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

base <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/homer2/'

run_homer <- function(){

	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')

	# limit it to unique SNPs in each method
	egenes_master <- egenes_master[pval_nominal < 0.01]
	egenes_master[,pair := paste0(variant_id,':',gene_id)]

	l <- egenes_master[method=='LAVA']
	g <- egenes_master[method=='global']

	repeats <- l[pair %in% g[,pair], gene_id]
	l <- l[!(gene_id %in% repeats)]
	g <- g[!(gene_id %in% repeats)]

	# now this includes only unique eGenes and same eGenes with different lead SNPs

	egenes <- data.table(rbind(l,g))
	egenes[,chr := sapply(variant_id, function(x) unname(unlist(strsplit(x,'_')))[1] )]
	egenes[,pos := sapply(variant_id, function(x) as.numeric(unname(unlist(strsplit(x,'_')))[2]) )]
	egenes[,start := pos-15]
	egenes[,stop := pos+15]
	egenes[,filler := NA]
	egenes[,strand := '+']

	for (t in tissues){

		outdir1 <- paste0(base, t, '/gvl')
		outdir2 <- paste0(base, t, '/lvg')
		system(sprintf('mkdir -p %s',outdir1))
		system(sprintf('mkdir -p %s',outdir2))

		global_bed <- egenes[method=='global' & tissue==t, .(chr,start,stop,variant_id,filler,strand)]
		global_bed <- global_bed[!duplicated(global_bed)]
		local_bed <- egenes[method=='LAVA' & tissue==t, .(chr,start,stop,variant_id,filler,strand)]
		local_bed <- local_bed[!duplicated(local_bed)]	

		out_local <- paste0(base, t, '/egenes_master.local.bed')
		out_global <- paste0(base, t, '/egenes_master.global.bed')

		local_sorted <- paste0(base, t, '/egenes_master.local.sorted.bed')
		global_sorted <- paste0(base, t, '/egenes_master.global.sorted.bed')

		write.table(global_bed, out_global, sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
		write.table(local_bed, out_local, sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

		system(sprintf('bedtools sort -i %s > %s',out_local,local_sorted))
		system(sprintf('bedtools sort -i %s > %s',out_global,global_sorted))

		system(sprintf('taskset -c 0-11 findMotifsGenome.pl %s hg38 %s -p 12 -bg %s &',global_sorted, outdir1, local_sorted))
		system(sprintf('taskset -c 0-11 findMotifsGenome.pl %s hg38 %s -p 12 -bg %s &',local_sorted, outdir2, global_sorted))

	}

}

# now analyze the outputs

analyze_homer <- function(tissue){

	dflist <- list()
	i <- 1
	
	system('echo "--------------------------------"')
	system(sprintf('echo %s',tissue))
	system('echo')
	system('echo "global v local:"')

	if(file.exists(paste0(base, tissue, '/gvl/knownResults.txt'))){
		
		gvl_known <- fread(paste0(base, tissue, '/gvl/knownResults.txt'), sep='\t', header=TRUE)
		gvl_known <- gvl_known[order(`q-value (Benjamini)`)]
		gvl_known[,target := 'Global']
		gvl_known[,background := 'Local']

		if(nrow(gvl_known[`q-value (Benjamini)` < 0.1]) == 0){
			system('echo "no motifs with FDR < 10% in global v local"')
		} else {
			print(gvl_known[`q-value (Benjamini)` < 0.1])
			dflist[[i]] <- gvl_known[`q-value (Benjamini)` < 0.1]
			i <- i + 1
		}

	} else {
		system('echo "no motifs with FDR < 10% in global v local"')
	}

	system('echo')
	system('echo "local v global:"')

	if(file.exists(paste0(base, tissue, '/lvg/knownResults.txt'))){
	
		lvg_known <- fread(paste0(base, tissue, '/lvg/knownResults.txt'), sep='\t', header=TRUE)
		lvg_known <- lvg_known[order(`q-value (Benjamini)`)]
		lvg_known[,target := 'Local']
		lvg_known[,background := 'Global']

		if(nrow(lvg_known[`q-value (Benjamini)` < 0.1]) == 0){
			system('echo "no motifs with FDR < 10% in local v global"')
		} else {
			print(lvg_known[`q-value (Benjamini)` < 0.1])
			dflist[[i]] <- lvg_known[`q-value (Benjamini)` < 0.1]
			i <- i + 1
		}

	} else {
		system('echo "no motifs with FDR < 10% in local v global"')
	}

	system('echo')

	results <- data.table(rbindlist(dflist))
	return(results)

}

#run_homer()
dtlist <- list()
i <- 1
for (t in tissues){
	print(t)
	results <- analyze_homer(t)
	if(nrow(results)==0){
		next
	}
	results[,tissue:=t]
	dtlist[[i]] <- results
	i <- i+1
}

results <- data.table(rbindlist(dtlist))
write.table(results,paste0(base,'homer-merged-results-fdr-10.txt'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)


# The true output of HOMER are "*.motif" files which contain the information necessary to identify future instance of motifs.  
# They are reported in the output directories from findMotifs.pl and findMotifsGenome.pl.  
# A typical motif file will look something like:

# >ASTTCCTCTT     1-ASTTCCTCTT    8.059752        -23791.535714   0       T:17311.0(44 ...
# 0.726   0.002   0.170   0.103
# 0.002   0.494   0.354   0.151
# 0.016   0.017   0.014   0.954
# 0.005   0.006   0.027   0.963
# 0.002   0.995   0.002   0.002
# 0.002   0.989   0.008   0.002
# 0.004   0.311   0.148   0.538
# 0.002   0.757   0.233   0.009
# 0.276   0.153   0.030   0.542
# 0.189   0.214   0.055   0.543

# The first row starts with a ">" followed by various information, and the other rows are the positions specific probabilities for each nucleotide (A/C/G/T).  
# The header row is actually TAB delimited, and contains the following information:

# ">" + Consensus sequence (not actually used for anything, can be blank) example: >ASTTCCTCTT
# Motif name (should be unique if several motifs are in the same file) example: 1-ASTTCCTCTT  or NFkB
# Log odds detection threshold, used to determine bound vs. unbound sites (mandatory) example: 8.059752
# log P-value of enrichment, example: -23791.535714
# 0 (A place holder for backward compatibility, used to describe "gapped" motifs in old version, turns out it wasn't very useful :)
# Occurence Information separated by commas, example: T:17311.0(44.36%),B:2181.5(5.80%),P:1e-10317
# 	T:#(%) - number of target sequences with motif, % of total of total targets
# 	B:#(%) - number of background sequences with motif, % of total background
# 	P:# - final enrichment p-value
# Motif statistics separated by commas, example: Tpos:100.7,Tstd:32.6,Bpos:100.1,Bstd:64.6,StrandBias:0.0,Multiplicity:1.13
# 	Tpos: average position of motif in target sequences (0 = start of sequences)
# 	Tstd: standard deviation of position in target sequences
# 	Bpos: average position of motif in background sequences (0 = start of sequences)
# 	Bstd: standard deviation of position in background sequences
# 	StrandBias: log ratio of + strand occurrences to - strand occurrences.
# 	Multiplicity: The averge number of occurrences per sequence in sequences with 1 or more binding site.
