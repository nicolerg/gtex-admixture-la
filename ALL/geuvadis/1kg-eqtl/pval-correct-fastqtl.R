# not currently needed

library(data.table)
library(qvalue)

# add p-value correction to results from fastQTL (permutation pass)

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
f <- args[2]
infile <- paste0(indir,'/',f)

if(grepl("AFR",f)){
	outfile <- paste0(indir, '/', 'AFR.geu.egenes.5fdr.txt')
} else if (grepl("EUR",f)){
	outfile <- paste0(indir, '/', 'EUR.geu.egenes.5fdr.txt')
}

egene <- fread(input=sprintf("zcat %s",infile),sep=' ',header=TRUE)

# BH correction
egene[,qval := qvalue(pval_beta)$qvalues]

egene <- egene[ qval < 0.05 ]

print(f)
print(nrow(egene))

write.table(egene, outfile, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
