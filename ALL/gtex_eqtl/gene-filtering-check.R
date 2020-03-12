# Spot check rejected genes to see if they pass the thresholds in Whole Blood

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

genelist <- c("ENSG00000257246",
"ENSG00000102174",
"ENSG00000166391",
"ENSG00000214062",
"ENSG00000269113",
"ENSG00000266644",
"ENSG00000254092",
"ENSG00000197421",
"ENSG00000196167",
"ENSG00000118972",
"ENSG00000218358",
"ENSG00000236226",
"ENSG00000204581",
"ENSG00000244265",
"ENSG00000164112",
"ENSG00000257246",
"ENSG00000272456",
"ENSG00000219395",
"ENSG00000211682",
"ENSG00000138669",
"ENSG00000225715",
"ENSG00000129455",
"ENSG00000257941",
"ENSG00000177133",
"ENSG00000273085")

s <- "Whole Blood"

samples <- annot$SAMPID[grepl(s, annot$SMTSD)]
nsub <- norm[,c(1,2,na.omit(match(samples,colnames(norm)))),with=FALSE]
csub <- counts[,c(1,2,na.omit(match(samples,colnames(counts)))),with=FALSE]
if(ncol(nsub) < 72 | ncol(csub) < 72){
	print("Not enough samples in the tissue")
} else {
	for (gene in genelist){
		nr <- as.numeric(nsub[grep(gene, nsub$gene_id),-c(1:2)])
		cr <- as.numeric(csub[grep(gene, nsub$gene_id),-c(1:2)])
		if(length(grep("TRUE", cr > 6)) < 20){
			print("Count threshold not met")
			print(cr[cr>6])
		} else if(length(grep("TRUE", nr > 0.1)) < 20){
			print("RPKM expression threshold not met")
			print(nr[nr>0.1])
		} else {
			print(gene)
		}
	}
}

