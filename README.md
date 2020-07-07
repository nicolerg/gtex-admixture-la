# Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx

This repository provides code and resources used to perform the analyses presented in the bioRxiv preprint [Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx](https://www.biorxiv.org/content/10.1101/836825v1) by Nicole Gay, *et al*. The author of all code is Nicole Gay (nicolerg@stanford.edu) unless otherwise noted. Please see the manuscript for context.   

Note that code is primarily provided to replicate analyses in the manuscript. However, when applicable, updated or more efficient versions of the code are also provided.  

"117AX" refers to a subset of 117 admixed invididuals in the GTEx v8 release whose GTEx IDs are provided in [metadata/gtex-admixed0.9.txt](metadata/gtex-admixed0.9.txt). These scripts could be adapted for any subset of individuals.  

**Please see a subdirectory for detailed documentation for each step:** 

### Subdirectories:  
- [**Local ancestry inference with RFMix:**](rfmix) Use RFMix v1.5.4 and 1000 Genomes reference panels to perform local ancestry inference on 117AX   
- [**eQTL calling:**](eqtl) Perform eQTL calling with both Global Ancestry Adjustment (GlobalAA) and Local Ancestry Adjustment (LocalAA)  
- [**Colocalization:**](colocalization) Perform colocalization between eQTL summary statistics (GlobalAA and LocalAA) and each of 114 GWAS (see [this preprint](https://www.biorxiv.org/content/10.1101/814350v1))
- [**Global ancestry inference with ADMIXTURE:**](admixture) Compare global ancestry fractions estimated from RFMix to ADMIXTURE global ancestry inference  
- [**Annotations:**](annotation) Annotate SNPs or genes with characteristics including variance in genotype explained by local ancestry, variance in gene expression explained by ancestry, RSIDs, LD, etc.  
- [**Plotting:**](plotting) Generate plots and supplementary tables included in the manuscript  
