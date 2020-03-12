# Split GTEx gene_reads and gene_tpm files into tissue-specific gene_reads and gene_tpm files based on sample annotation
# Filter to subjects that also have genotype information

library(data.table)

readcounts <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/rna_seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
readnorm <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/rna_seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
annotfile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
genotlist <- "/users/nicolerg/gtex-admix/metadata/gtex-genot.txt" # header from VCF

counts <- fread(input = sprintf("zcat %s", readcounts), skip=2, header=TRUE, sep='\t', stringsAsFactors=FALSE)
norm <- fread(input = sprintf("zcat %s", readnorm), skip=2, header=TRUE, sep='\t', stringsAsFactors=FALSE)
annot <- read.delim(annotfile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
genot <- read.delim(genotlist, sep='\t', header=TRUE, stringsAsFactors=FALSE)
genot <- colnames(genot)
genot <- gsub(".","-",genot,fixed=TRUE)
genot <- genot[-c(1:9)]

# subset annot file by samples that also have genotype information
# SAMPID: GTEX-1117F-0003-SM-6WBT7
fun <- function(a){paste(unlist(strsplit(a,"-"))[1:2],collapse='-')}
annot$SUBJID <- sapply(annot$SAMPID,fun)

for (id in unique(annot$SUBJID)){
	if (!as.character(id)%in%as.character(genot)) {
		annot <- annot[!(annot$SUBJID==id),]
	}
}

sys <- unique(annot$SMTSD)
sys <- sys[sys!=""]

outdir <- "/mnt/lab_data/montgomery/nicolerg/gtex-sampleannot"

for (s in sys){
	samples <- annot$SAMPID[grepl(s, annot$SMTSD)]

	nsub <- norm[,c(1,2,na.omit(match(samples,colnames(norm)))),with=FALSE]
	row1 <- rep(NA, ncol(nsub))
	row1[1] <- "#1.2"
	row2 <- rep(NA, ncol(nsub))
	row2[1] <- nrow(nsub)
	row2[2] <- ncol(nsub) - 2
	row1 <- t(data.table(row1))
	row2 <- t(data.table(row2))
	row3 <- t(data.table(colnames(nsub)))
	colnames(row1) <- colnames(nsub)
	colnames(row2) <- colnames(nsub)
	colnames(row3) <- colnames(nsub)
	nsub <- rbind(row1, row2, row3, nsub)

	csub <- counts[,c(1,2,na.omit(match(samples,colnames(counts)))),with=FALSE]
	row1 <- rep(NA, ncol(csub))
	row1[1] <- "#1.2"
	row2 <- rep(NA, ncol(csub))
	row2[1] <- nrow(csub)
	row2[2] <- ncol(csub) - 2
	row1 <- t(data.table(row1))
	row2 <- t(data.table(row2))
	row3 <- t(data.table(colnames(csub)))
	colnames(row1) <- colnames(csub)
	colnames(row2) <- colnames(csub)
	colnames(row3) <- colnames(csub)
	csub <- rbind(row1, row2, row3, csub)

	s <- gsub(" - ","_",s,fixed=TRUE)
	s <- gsub(" (","_",s,fixed=TRUE)
	s <- gsub(")","_",s,fixed=TRUE)
	s <- gsub(" ","_",s,fixed=TRUE)
	write.table(nsub,paste(outdir,"/",s,"_gene_tpm.gct",sep=''),sep='\t',quote=FALSE,na='',row.names=FALSE,col.names=FALSE)
	write.table(csub,paste(outdir,"/",s,"_gene_reads.gct",sep=''),sep='\t',quote=FALSE,na='',row.names=FALSE,col.names=FALSE)
}