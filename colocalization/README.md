# Colocalization

## Download GWAS summary statistics 
GWAS summary statistics were downloaded from two sources: 
  - 28 GWAS from the PAGE Consortium, presented in ["Genetic analyses of diverse populations improves discovery for complex traits"](https://www-nature-com.stanford.idm.oclc.org/articles/s41586-019-1310-4), available through the [GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics)
  - 114 GWAS from ["Widespread dose-dependent effects of RNA expression and splicing on complex diseases and traits"](https://www.biorxiv.org/content/10.1101/814350v1) (available [here](https://zenodo.org/record/3629742#.Xn-eedNKhQJ))

## Format the GWAS summary statistics 
[`format_gtex_gwas.sh`](format_coloc.sh) was used to convert the 114 imputed GWAS summary statistics to the format required by the coloc wrapper pipeline.  

[`format_page_gwas.sh`](format_page_gwas.sh) was used to covert the PAGE GWAS summary statistics to the format required by the coloc wrapper pipeline.  

## Perform colocalization with wrapper pipeline 

### COLOC
Follow [`coloc_pipeline.sh`](coloc_pipeline.sh) to do the following:  
  1. Calculate SNP allele frequencies for all tested SNPs based on 117AX genotypes (`snp_to_effect_af.tsv.gz`)  
  2. Parse allpairs files to include only tests for "same eGene, different lead SNP" genes at a nominal p-value of 1e-04. This could be optimized by making these gene lists tissue-specific. More tests are currently performed than necessary. See [`parse_allpairs_for_coloc.py`](parse_allpairs_for_coloc.py)  
  3. Sort, `bgzip`, and `tabix`-index filterd allpairs files  
  4. Make config files for the colocalization wrapper pipeline (see [`gwas_char.txt`](config/gwas_char.txt), [`gwas_experiments.json`](config/gwas_experiments.json), [`gtex_cc.json`](gtex_cc.json), [`format_config.py`](format_config.py); all config files used for COLOC are available [here](config/coloc))  
  5. Make `n_gwas.txt` and `n_eqtl.txt` sample size files manually (or include sample sizes in the config file)  
  6. Run COLOC with the colocalization wrapper pipeline  
  7. Concatenate results  

Here is the code used to generate colocalization posterior probabilities with COLOC:
```r
suppressMessages(require(coloc))

# Get input and output locations as command-line arguments
args = commandArgs(TRUE)
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

### FINEMAP  
We performed FINEMAP colocalization on the subset of loci for which COLOC provided evidence of colocalization (PP4 > 0.5).  
1. Adjust the COLOC outputs to ensure that the same lead SNP is used to test both eQTL files per tissue (see [`adjust_coloc_sites.py`](adjust_coloc_sites.py))  
2. Filter sites by COLOC PP4, i.e. `signif_coloc_loci.tests.txt` (see [`filter_loci_for_finemap.R`](filter_loci_for_finemap.R))  
3. Run FINEMAP (v1.1) with the colocalization wrapper pipeline (see [`run_finemap_sig.sh`](run_finemap_sig.sh) and [config file](config/finemap/finemap_1000G_sig.json))   

Within the wrapper pipeline, FINEMAP is called with the following parameters: 
```bash
finemap --sss --in-files ${infile} --n-causal-max 1 --n-iterations 1000000 --n-convergence 1000
```
For each locus, the colocalization posterior probability (CLPP) score is calculated from the marginal Posterior Inclusion Probabilities [as previously described](https://www-ncbi-nlm-nih-gov.stanford.idm.oclc.org/pmc/articles/PMC5142122/). 

**NOTE:** To run the colocalization wrapper pipeline on SCG, use the `scg` branch of the `brain_gwas` repository and load the following modules first IN THIS ORDER: 
```bash 
module load r/3.6 # R v3.6
module load miniconda/2 # python2
module load plink/1.90b6.13 # PLINK v1.90b
module load tabix # tabix
```

You will also need to install [`FINEMAP v1.1`](http://www.christianbenner.com/), soft-link the executable to a file named `finemap`, and add the directory with the `finemap` soft link to your PATH. 

Contact Mike Gloudemans about access to the colocalization wrapper pipeline: mgloud@stanford.edu

## Merge colocalization results from both methods  
After concatenating the colocalization results from each method, merge the results. [`format_coloc.R`](format_coloc.R) provides the code for how this was done specifically for our project (the results are saved in `master_coloc-1e-04-20200418.RData`). 

## Create colocalization signal plots  
[`locus_plots.R`](locus_plots.R) generates tables used to make colocalization signal plots for loci in which one ancestry adjustment method has a stronger colocalization than the other. These plots are shown in Figure S4 of [our original manuscript](https://www.biorxiv.org/content/10.1101/836825v1) and Figure S5 of our revised manuscript. 
