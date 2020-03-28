# Colocalization

## Download GWAS summary statistics 
GWAS summary statistics were downloaded from two sources: 
  - 28 GWAS from the PAGE Consortium, presented in ["Genetic analyses of diverse populations improves discovery for complex traits"](https://www-nature-com.stanford.idm.oclc.org/articles/s41586-019-1310-4), available through the [GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics)
  - 114 GWAS from ["Widespread dose-dependent effects of RNA expression and splicing on complex diseases and traits"](https://www.biorxiv.org/content/10.1101/814350v1) (available [here](https://zenodo.org/record/3629742#.Xn-eedNKhQJ))

## Format the GWAS summary statistics 
[`format_gtex_gwas.sh`](format_coloc.sh) was used to convert the 114 "Widespread..." GWAS summary statistics to the format required by the coloc wrapper pipeline.  

[`format_page_gwas.sh`](format_page_gwas.sh) was used to covert the PAGE GWAS summary statistics to the format required by the coloc wrapper pipeline.  

## Run the colocalization pipeline (COLOC and FINEMAP)  
Follow [`coloc_pipeline.sh`](coloc_pipeline.sh) to do the following:  
  1. Calculate SNP allele frequencies for all tested SNPs based on 117AX genotypes (`snp_to_effect_af.tsv.gz`) 
  2. [`parse_allpairs_for_coloc.py`](parse_allpairs_for_coloc.py): parse allpairs files to include only tests for "same eGene, different lead SNP" genes at a nominal p-value of 1e-04. (This could be optimized by making these gene lists tissue-specific. More tests are currently performed than necessary.)  
  3. Sort, `bgzip`, and `tabix`-index filterd allpairs files  
  4. Make config file for the coloc wrapper pipeline (see [`gwas_char.txt`](config/gwas_char.txt) and [`config.json`](config/config.json); all config files used are available [here](config))    
  5. Make `n_gwas.txt` and `n_eqtl.txt` sample size files manually (or include sample sizes in the config file)
  6. Run COLOC and FINEMAP with the colocalization wrapper pipeline 
  7. Merge results  

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
