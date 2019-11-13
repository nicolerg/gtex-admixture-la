# filter gene expression down to 117AX

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
expr <- args[1] # GTEx normalized expression file 
outdir <- args[2]

keep = fread("gtex-admixed0.9.txt", sep='\t', header=F) # see "metadata" in repository
expression = fread(cmd=sprintf('zcat %s',expr),sep='\t',header=T)

# filter columns 
keep_cols <- c('#chr','start','end','gene_id',keep[,V1])
expression = expression[,c(colnames(expression)[colnames(expression) %in% keep_cols]),with=FALSE]

tissue <- unlist(strsplit(basename(expr), '[.]'))[1]
write.table(expression, sprintf('%s/%s.v8.normalized_expression.admixed_subset.bed',outdir,tissue), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
