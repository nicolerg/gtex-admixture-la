# filter gene expression down to 117AX

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
expr <- args[1] # GTEx normalized expression file 
subj <- args[2]
outdir <- args[3]

keep = fread(subj, sep='\t', header=F) # see root of repository: "gtex-admixed0.9.txt"
expression = fread(cmd=sprintf('zcat %s',expr),sep='\t',header=T)

# filter columns 
keep_cols <- c('#chr','start','end','gene_id',keep[,V1])
expression = expression[,c(colnames(expression)[colnames(expression) %in% keep_cols]),with=FALSE]

tissue <- unlist(strsplit(basename(expr), '[.]'))[1]
write.table(expression, sprintf('%s/%s.v8.normalized_expression.admixed_subset.bed',outdir,tissue), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
