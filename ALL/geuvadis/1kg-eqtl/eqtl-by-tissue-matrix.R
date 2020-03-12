library(data.table)

# genes <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/gtex-genelist.txt"
# g <- fread(genes, sep='\t', header=FALSE)
# colnames(g) <- "gene_id"
# genes <- g[,gene_id]

# path <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/"
# file.names <- dir(path, pattern ="v8.egenes.txt.gz")

# master <- data.table(gene_id=genes)

# for(f in file.names){
# 	input <- paste0(path,f)
# 	egenes <- fread(input=sprintf("zcat %s",input), header=TRUE, sep="\t")
# 	system <- gsub(path,"",f)
# 	system <- gsub(".v8.egenes.txt.gz","",system,fixed=TRUE)
# 	print(system)

# 	eqtlgene <- egenes[ pval_nominal <= pval_nominal_threshold, gene_id ]
# 	sigbin <- rep(NA, length(genes))
# 	for(i in 1:length(genes)){
# 		g <- genes[i]
# 		sigbin[i] <- as.integer(g %in% eqtlgene)
# 	}
# 	print(sigbin)
# 	master[,(system):=sigbin]
# }

# write.table(master, "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/gtex-eqtl-genes-bytissue.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

# instead, look to see if those genes are even included in the normalized expression file for the tissue

genes <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/eqtl.curated.geneset.chr22.tsv"
g <- fread(genes, sep='\t', header=TRUE)
genes <- g[,pid]
genes <- gsub("[.].*","",genes)
print(genes)

path <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/"
file.names <- dir(path, pattern ="v8.normalized_expression.bed.gz")
file.names <- file.names[!grepl(".gz.tbi",file.names)]

master <- data.table(gene_id=genes)

for(f in file.names){
	input <- paste0(path,f)
	expr <- fread(input=sprintf("zcat %s",input), header=TRUE, sep="\t")
	
	expr[,gene_id:=gsub("[.].*","",gene_id)]
	print(expr[1:5,gene_id])

	system <- gsub(path,"",f)
	system <- gsub(".v8.normalized_expression.bed.gz","",system,fixed=TRUE)
	print(system)

	sigbin <- rep(NA, length(genes))
	i <- 1
	for (g in genes){
		sigbin[i] <- as.integer(g %in% expr[,gene_id])
		i <- i + 1
	}
	print(sigbin)
	master[,(system):=sigbin]
}

write.table(master, "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/gtex-genes-expressed-bytissue-chr22.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

master <- fread("/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/gtex-genes-expressed-bytissue-chr22.tsv", sep='\t', header=TRUE)
genes <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/eqtl.curated.geneset.chr22.tsv"
g <- fread(genes, sep='\t', header=TRUE)
genes <- g[,pid]
genes <- gsub("[.].*","",genes)

m <- na.omit(match(genes, master[,gene_id]))
matfilt <- master[m]
matfilt[,gene_id:=NULL]
df <- data.frame(sort(colSums(matfilt)/length(genes), decreasing=TRUE))
colnames(df) <- "frac_expr"
df$rank <- 1:nrow(df)
df
