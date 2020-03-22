# Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx

This repository provides code and resources used to perform the analyses presented in the bioRxiv preprint [Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx](https://www.biorxiv.org/content/10.1101/836825v1) by Nicole Gay, et al. The author of all code is Nicole Gay (nicolerg@stanford.edu) unless otherwise noted. Please see the manuscript for context.   

Note that code is primarily provided to replicate analyses in the manuscript. However, when applicable, updated or more efficient versions of the code are also provided.  

"117AX" refers to a subset of 117 admixed invididuals in the GTEx v8 release whose GTEx IDs are provided in [`metadata/gtex-admix0.9.txt`](metadata/gtex-admix0.9.txt). These scripts could be adapted for any subset of individuals.  

## Table of contents  
- **Perform local ancestry inference:** Use RFMix v1.5.4 and 1000 Genomes reference panels to perform local ancestry inference on 117AX   
- **eQTL calling:** Perform eQTL calling with both Global Ancestry Adjustment (GlobalAA) and Local Ancestry Adjustment (LocalAA)  
- **Colocalization:** <Add something here>

## eQTL calling

### Prepare VCFs 

The SHAPEIT2-phased, MAF 0.01-filtered GTEx v8 WGS VCF was converted to hg19 using LiftOver for compatibility with the GRCh37 HapMap genetic map and hg19 1000 Genomes ([`liftOver-hg38Tohg19.sh`](rfmix/liftOver-hg38Tohg19.sh)). This VCF was then split into chromosome-level VCFs, each of which was then filtered to include only 117AX samples and MAC > 10 using `vcftools`:
```
vcfin=genotypes-hg19
vcfdir=admixed/geno

for chr in {1..22}; do
	vcf_unfilt=${vcfin}/gtex.MAF01.phased.chr${chr}.vcf.gz
	vcf_out=${vcfdir}/gtex.admixed.MAC10.phased.chr${chr}
	vcftools \
		--gzvcf ${vcf_unfilt} \
		--keep ${admix_ids} \
		--recode-INFO-all \
		--mac 10 \
		--max-alleles 2 \
		--min-alleles 2 \
		--out ${vcf_out} \
		--recode 
done

gzip ${vcfdir}/gtex.admixed.MAC10.phased.chr*
```

### Prepare covariates

[`prepare_covariates.sh`](eqtl/prepare_covariates.sh) concatenates covariates provided by GTEx and new PEER factors. It requires paths to two GTEx v8 files:  
  - `${tissue}.v8.normalized_expression.bed.gz` (hg38 coordinates now)   
  - `${tissue}.v8.covariates.txt`  

It also calls 3 other scripts:  
  - [`filter_expression_admixed.R`](eqtl/filter_expression_admixed.R)
  - [`concat_cov.R`](eqtl/concat_cov.R)
  - [`run_PEER.R`](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R) (not in this repository)

It has a few outputs:  
  - `${tissue}.v8.normalized_expression.admixed_subset.bed`: expression in 117AX (however many have samples in this tissue)
  - `other_covariates.txt`: necessary covariates extracted from the GTEx v8 covariates file  
  - `peer` folder with 3 outputs, including `Whole_Blood.PEER_covariates.txt`, which are the PEER covariates calculated using only 117AX samples
  - `${tissue}.all_covariates.txt`: combined covariates. **This file will be used for eQTL calling.**

### Run eQTL calling with both LocalAA and GlobalAA 

[`batch_eqtl_localaa_globalaa.sh`](eqtl/batch_eqtl_localaa_globalaa.sh) splits the expression file for a tissue into chunks of 50 genes (line 49) and runs [`eqtl_localaa_globalaa.R`](eqtl/eqtl_localaa_globalaa.R) on each chunk of genes per chromosome.  

[`eqtl_localaa_globalaa.R`](eqtl/eqtl_localaa_globalaa.R) requires several inputs (see arguments):
  - `chr`: Chromosome number
  - `tissue`: Tissue specified with GTEx file prefix, e.g. `Muscle_Skeletal` 
  - `exprfile`: Path to GTEx v8 expression file for this tissue 
  - `globalcov`: Path to `${tissue}.all_covariates.txt` from the [previous step](#prepare-covariates)
  - `geno`: Path to filtered VCF for this chromosome (from [this step](#prepare-vcfs))
  - `localcov`: Path to local ancestry covariates for this chromosome from [this step](#make-master-local-ancestry-files)
  - `out`: Output file. Results from each test are added to this file as they are generated. This file is analogous to GTEx's "allpairs" files  

[`batch_eqtl_localaa_globalaa.sh`](eqtl/batch_eqtl_localaa_globalaa.sh) was written to parallelize eQTL calling on a cluster without a job submission system. As written, it restricts processes to 9 cores. If you have access to a HPC with a job submission system like `SGE` or `SLURM`, I highly recommend parallelizing this process as much as possible. See scripts in [eqtl/slurm](eqtl/slurm) for an example of how to do this with `SLURM`.  

[`eqtl_localaa_globalaa.R`](eqtl/eqtl_localaa_globalaa.R) generates a progress bar for its corresponding chunk for your sanity. Once all chunks are done running, concatenate the results for each tissue (see the end of [`batch_eqtl_localaa_globalaa.sh`](eqtl/batch_eqtl_localaa_globalaa.sh)).  

Now let's assume that, for each tissue, you have a master allpairs file called `${tissue}-LocalAA-GlobalAA-allpairs-merged.tsv.gz`. These files are very large. Some of the downstream steps are facilitated by having a filtered allpairs file that only includes tests where either the GlobalAA or LocalAA nominal p-value is less than 0.05. Generate this filtered allpairs file as follows:  
```
zcat ${tissue}-LocalAA-GlobalAA-allpairs-merged.tsv.gz | \
  awk -F "\t" '{ if(($6 < 0.05) || ($9 < 0.05)) { print } }' \
  > ${tissue}-LocalAA-GlobalAA-allpairs-merged-filt.tsv
gzip ${tissue}-LocalAA-GlobalAA-allpairs-merged-filt.tsv
```

### Get tied lead SNPs for each gene 

[`extract-egenes.py`](eqtl/extract-egenes.py) identifies the most significant SNP(s) for each tested gene. If multiple SNPs have the same smallest p-value, all tied lead SNPs are reported. No significance threshold is applied; the lead SNP(s) is reported for **all** tested genes.  

The input should be the **unfiltered** allpairs file for a given tissue. There are two outputs:
  - `gtex.admix.global.egenes.tied.txt.gz`: Lead SNPs according to GlobalAA  
  - `gtex.admix.lava.egenes.tied.txt.gz`: Lead SNPs according to LocalAA  

### Generate some eGene sets that are repeatedly used in downstream analyses





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
