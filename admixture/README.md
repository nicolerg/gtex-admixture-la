# Global ancestry inference with ADMIXTURE

Correlate GTEx genotype PCs with RFMix [global ancestry estimates](https://github.com/nicolerg/gtex-admixture-la/tree/master/rfmix#calculate-global-ancestry-proportions):
```{r}
library(data.table)
library(ggcorrplot)

pcs <- fread('GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',sep='\t',header=T)
global_ai <- fread('gtex_ai_global_all.txt',sep='\t',header=T) # global ancestry proportions from RFMix
pcs[,SUBJID := sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
m <- merge(pcs, global_ai, by.x='SUBJID', by.y='ID')
m[,UNK := NULL]
m <- data.frame(m)
rownames(m) <- m$SUBJID
m$SUBJID <- NULL
m$FID <- NULL
m$IID <- NULL
corr <- cor(m)
ggcorrplot(corr[1:10, 21:23],lab=T)
```

Use `bcftools` to merge GTEx and 1000 Genomes VCFs (hg38). Generate one merged, MAF-filtered VCF per chromosome. [Snakefile](Snakefile) provides a way to parallellize this using a Snakemake pipeline. These VCFs currently include all samples from both consortia.  

ADMIXTURE was run with several different combinations of reference panels:
1. Same reference panels as used in RFMix 
2. Same superpopulations as used in RFMix but only 100 reference panels from each superpopulation (equal numbers from each subpopulation) 
3. Same as (2) but with the addition of SAS 

ADMIXTURE was also run on 
