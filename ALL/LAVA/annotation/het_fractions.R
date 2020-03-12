library(data.table)
library(reshape2)

outdir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/het_fraction'
infile <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/het_fraction/diff-lead-pairs-by-tissue-1e-04.txt'
edir <- '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices'
admixed <- fread('/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt',sep='\t',header=F)
admixed <- admixed[,V1]

pairs <- fread(infile, sep='\t', header=T)
pairs <- pairs[,.(variant_id_local,variant_id_global,gene_id,tissue)]
pairs <- data.table(melt(pairs, measure.vars=c('variant_id_global','variant_id_local'), variable.name='method',value.name='variant_id'))

pairs[,chrom := sapply(variant_id, function(x) unname(unlist(strsplit(x, '_')))[1])]

# initialize output file
hetfrac <- file('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/het_fraction/diff-lead-pairs-by-tissue-1e04-het-frac.tsv', 'w')
header=paste('variant_id','gene','tissue','n_subj','het_frac',sep='\t')
cat(paste0(header,'\n'), file=hetfrac)

# iterate through chromosomes and output fraction of het chromosomes for each SNP (load each ancestry file only once)
for (chr in rev(1:22)){

	print(chr)

	# filter pairs down to this chromosome
	pairs_chr <- pairs[chrom == paste0('chr',chr)]

	# load the local ancestry file
	anc <- fread(cmd=sprintf('zcat /mnt/lab_data/montgomery/nicolerg/admixed-ancestry/chr%s.hg19.maf0.localcov.tsv.gz',chr), sep='\t', header=T)

	for (t in c('Muscle_Skeletal','Lung','Artery_Tibial','Nerve_Tibial','Adipose_Subcutaneous','Skin_Not_Sun_Exposed_Suprapubic')){

		sub <- pairs_chr[tissue==t]

		if(!paste0('chr',chr)%in%sub[,chrom]){
			next
		}

		# need to get the subjects in that tissue 
		expr <- fread(cmd=sprintf('zcat %s/%s.v8.normalized_expression.bed.gz',edir,t), sep='\t', header=T, nrows=1)
		subj_in_tissue <- colnames(expr)[5:ncol(expr)] 

		# filter down to subjects in this tissue
		anc_t <- anc[ (SUBJID %in% subj_in_tissue) & (SUBJID %in% admixed) ]

		for (snp in sub[,variant_id]){

			# filter down to this SNP
			snp_anc <- anc_t[SNP_ID == snp]
			
			# remove subjects where UNK is not 0
			snp_anc <- snp_anc[UNK==0]

			# get fraction of het subjects
			total_subj <- length(subj_in_tissue)
			valid_subj <- nrow(snp_anc)
			snp_anc[,HET := ifelse(EUR==1 | AFR==1 | ASN==1, 1, 0)]
			het_subj <- sum(snp_anc[,HET])
			frac_het <- het_subj/valid_subj

			# print to file
			gene <- sub[variant_id == snp, gene_id]

			string=paste(snp, gene, t, valid_subj, frac_het,
				sep='\t')
			writeLines(string)
			cat(paste0(string,'\n'), file=hetfrac, append=TRUE)

		}
	}
}

close(hetfrac)
