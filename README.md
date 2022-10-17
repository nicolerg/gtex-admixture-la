# Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx

This repository provides code and resources used to perform the analyses presented in the bioRxiv preprint [Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02113-0) by Nicole Gay, *et al*. The author of all code is Nicole Gay (nicolerg@stanford.edu) unless otherwise noted. Please see the manuscript for context.   

Note that code is primarily provided to replicate analyses in the manuscript. However, when applicable, updated or more efficient versions of the code are also provided.  

"117AX" refers to a subset of 117 admixed invididuals in the GTEx v8 release whose GTEx IDs are provided in [metadata/gtex-admixed0.9.txt](metadata/gtex-admixed0.9.txt). These scripts could be adapted for any subset of individuals.  

**Please see a subdirectory for detailed documentation for each step:** 

### Subdirectories:  
- [**Local ancestry inference with RFMix:**](rfmix) Use RFMix v1.5.4 and 1000 Genomes reference panels to perform local ancestry inference on 117AX   
- [**eQTL calling:**](eqtl) Perform eQTL calling with both Global Ancestry Adjustment (GlobalAA) and Local Ancestry Adjustment (LocalAA)  
- [**Colocalization:**](colocalization) Perform colocalization between eQTL summary statistics (GlobalAA and LocalAA) and each of 114 GWAS (see [this preprint](https://www.biorxiv.org/content/10.1101/814350v2))
- [**Global ancestry inference with ADMIXTURE:**](admixture) Compare global ancestry fractions estimated from RFMix to ADMIXTURE global ancestry inference  
- [**Annotations:**](annotation) Annotate SNPs or genes with characteristics including variance in genotype explained by local ancestry, variance in gene expression explained by ancestry, RSIDs, LD, etc.  
- [**Plotting:**](plotting) Generate plots and supplementary tables included in the manuscript  

## Citing 

**Manuscript:** Gay, N.R., Gloudemans, M., Antonio, M.L. *et al.* Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx. *Genome Biol* 21, 233 (2020). https://doi.org/10.1186/s13059-020-02113-0 

**Code:** Nicole R. Gay. (2020). Source code: Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.3924788

**Extended data:** Nicole R. Gay, Michael Gloudemans, Margaret L. Antonio, Brunilda Balliu, YoSon Park, Alicia R. Martin, Shaila Musharoff, Abhiram Rao, François Aguet, Alvaro Barbeira, Rodrigo Bonazzola, Farhad Hormozdiari, GTEx Consortium, Kristin G. Ardlie, Christopher D. Brown, Hae Kyung Im, Tuuli Lappalainen, Xiaoquan Wen, & Stephen B. Montgomery. (2020). Extended data: Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3926871

