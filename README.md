# Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx

This repository provides code and resources used to perform the analyses presented in the bioRxiv preprint [Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx](https://www.biorxiv.org/content/10.1101/836825v1) by Nicole Gay, et al. The author of all code is Nicole Gay (nicolerg@stanford.edu) unless otherwise noted. Please see the manuscript for context.   

Note that code is primarily provided to replicate analyses in the manuscript. However, when applicable, updated or more efficient versions of the code are also provided.  

"117AX" refers to a subset of 117 admixed invididuals in the GTEx v8 release whose GTEx IDs are provided in [`metadata/gtex-admixed0.9.txt`](metadata/gtex-admixed0.9.txt). These scripts could be adapted for any subset of individuals.  

## Table of contents  
- **Perform local ancestry inference:** Use RFMix v1.5.4 and 1000 Genomes reference panels to perform local ancestry inference on 117AX   
- **eQTL calling:** Perform eQTL calling with both Global Ancestry Adjustment (GlobalAA) and Local Ancestry Adjustment (LocalAA)  
- **Colocalization:** <Add something here>






---------------------------------

#TODO

## Colocalization (colocalization subdirectory)

< need to add get_uniq_egenes.R >

`coloc_pipeline.sh`
calls:
  - `parse_allpairs_for_coloc.py`: extract subset of eQTL allpairs needed for colocalization tests 
  - `compile_coloc.R`: merge coloc results 

Contact Mike Gloudemans about access to the colocalization wrapper pipeline: mgloud@stanford.edu

Here is the code used to generate colocalization posterior probabilities with COLOC:
```r
suppressMessages(require(coloc))

# Get input and output locations as command-line arguments
args <- commandArgs(TRUE)
infile = args[1]
N_gwas = as.numeric(args[2]) # sample size of GWAS
s_gwas = as.numeric(args[3])
type_gwas = args[4]
N_eqtl = as.numeric(args[5]) # sample size of eQTL study

# Run coloc
data = read.table(infile, header=TRUE, sep=',', fill=TRUE)

# remove incomplete cases
data = data[complete.cases(data),]

# GWAS data
dataset1 = list(pvalues=data$pvalue_gwas, MAF=data$effect_af_gwas, s=s_gwas, N=N_gwas, type=type_gwas)
# eQTL data
dataset2 = list(pvalues=data$pvalue_eqtl, MAF=data$effect_af_eqtl, N=N_eqtl, type="quant")

sink(file="/dev/null")
results = coloc.abf(dataset1, dataset2)
sink()

#summary is a vector giving the number of SNPs analysed,
#and the posterior probabilities of H0 (no causal variant),
#H1 (causal variant for trait 1 only),
#H2 (causal variant for trait 2 only),
#H3 (two distinct causal variants) and
#H4 (one common causal variant)

h0 = results$summary[2]
h1 = results$summary[3]
h2 = results$summary[4]
h3 = results$summary[5]
h4 = results$summary[6]

cat(h0, h1, h2, h3, h4, sep='\t')
```
`format_coloc.R`: clean up the colocalization results

`locus_plots.R`: generate tables used to make colocalization signal plots for loci in which one ancestry adjustment method has a stronger colocalization than the other; used to make Figure S4.
