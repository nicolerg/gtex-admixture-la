# Global ancestry inference with ADMIXTURE

Use `bcftools` to merge GTEx and 1000 Genomes VCFs (hg38). Generate one merged, MAF-filtered VCF per chromosome. [Snakefile](Snakefile) provides a way to parallellize this using a Snakemake pipeline. These VCFs currently include all samples from both consortia.  

Download ADMIXTURE v1.3.0: http://software.genetics.ucla.edu/admixture/download.html 

Then follow [`run_admixture.sh`](run_admixture.sh) to prepare variant files and run ADMIXTURE for a given set of reference panels (1000 Genomes samples). As written, this script runs ADMIXTURE using the same reference panels used to estimate local ancestry with RFMix, but it was run with several different combinations of reference panels:
1. Same reference panels as used in RFMix 
2. Same superpopulations as used in RFMix but only 100 reference panels from each superpopulation (equal numbers from each subpopulation) 
3. Same as (2) but with the addition of SAS 

[`run_admixture.sh`](run_admixture.sh) calls two other scripts:
  - [select_reference_populations.R](select_reference_populations.R) is used to select the 1000 Genomes reference panels. Modify this script as necessary (for example, to select equal numbers of reference panels from each superpopulation)  
  - [order_pop_file.R](order_pop_file.R) makes sure that the order of samples in the `fam` and `pop` files are the same 
  
 [`plot_admixture.R`](plot_admixture.R) provides an example of how to generate plots to compare ADMIXTURE global ancestry estimates, RFMix global ancestry estimates, gPCs, and GTEx self-reported ancestry. 
