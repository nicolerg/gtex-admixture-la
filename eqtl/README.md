# eQTL calling

## Prepare VCFs 

The SHAPEIT2-phased, MAF 0.01-filtered GTEx v8 WGS VCF was converted to hg19 using LiftOver for compatibility with the GRCh37 HapMap genetic map and hg19 1000 Genomes ([`liftOver-hg38Tohg19.sh`](../rfmix/liftOver-hg38Tohg19.sh)). This VCF was then split into chromosome-level VCFs, each of which was then filtered to include only 117AX samples and MAC > 10 using `vcftools`:
```bash
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

## Prepare covariates

[`prepare_covariates.sh`](prepare_covariates.sh) concatenates covariates provided by GTEx and new PEER factors. It requires paths to two GTEx v8 files:  
  - `${tissue}.v8.normalized_expression.bed.gz` (hg38 coordinates now)   
  - `${tissue}.v8.covariates.txt`  

It also calls 3 other scripts:  
  - [`filter_expression_admixed.R`](filter_expression_admixed.R)
  - [`concat_cov.R`](concat_cov.R)
  - [`run_PEER.R`](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R) (not in this repository)

It has a few outputs:  
  - `${tissue}.v8.normalized_expression.admixed_subset.bed`: expression in 117AX (however many have samples in this tissue)
  - `other_covariates.txt`: necessary covariates extracted from the GTEx v8 covariates file  
  - `peer` folder with 3 outputs, including `${tissue}.PEER_covariates.txt`, which are the PEER covariates calculated using only 117AX samples
  - `${tissue}.all_covariates.txt`: combined covariates. **This file will be used for eQTL calling.**

## Run eQTL calling with both LocalAA and GlobalAA 

[`batch_eqtl_localaa_globalaa.sh`](batch_eqtl_localaa_globalaa.sh) splits the expression file for a tissue into chunks of 50 genes (line 49) and runs [`eqtl_localaa_globalaa.R`](eqtl_localaa_globalaa.R) on each chunk of genes per chromosome.  

[`eqtl_localaa_globalaa.R`](eqtl_localaa_globalaa.R) requires several inputs (see arguments):
  - `chr`: Chromosome number
  - `tissue`: Tissue specified with GTEx file prefix, e.g. `Muscle_Skeletal` 
  - `exprfile`: Path to GTEx v8 expression file for this tissue 
  - `globalcov`: Path to `${tissue}.all_covariates.txt` from the [previous step](#prepare-covariates)
  - `geno`: Path to filtered VCF for this chromosome (from [this step](#prepare-vcfs))
  - `localcov`: Path to local ancestry covariates for this chromosome from [this step](https://github.com/nicolerg/gtex-admixture-la/tree/master/rfmix#make-master-local-ancestry-files)
  - `out`: Output file. Results from each test are added to this file as they are generated. This file is analogous to GTEx's "allpairs" files  

[`batch_eqtl_localaa_globalaa.sh`](batch_eqtl_localaa_globalaa.sh) was written to parallelize eQTL calling on a cluster without a job submission system. As written, it restricts processes to 9 cores. If you have access to a HPC with a job submission system like `SGE` or `SLURM`, I highly recommend parallelizing this process as much as possible. See scripts in [eqtl/slurm](slurm) for an example of how to do this with `SLURM`.  

[`eqtl_localaa_globalaa.R`](eqtl_localaa_globalaa.R) generates a progress bar for its corresponding chunk for your sanity. Once all chunks are done running, concatenate the results for each tissue (see the end of [`batch_eqtl_localaa_globalaa.sh`](batch_eqtl_localaa_globalaa.sh)).  

Now let's assume that, for each tissue, you have a master allpairs file called `${tissue}-LocalAA-GlobalAA-allpairs-merged.tsv.gz`. These files are very large. Some of the downstream steps are facilitated by having a filtered allpairs file that only includes tests where either the GlobalAA or LocalAA nominal p-value is less than 0.05. Generate this filtered allpairs file as follows:  
```bash
zcat zcat ${tissue}-LocalAA-GlobalAA-allpairs-merged.tsv.gz | head -1 > ${tissue}-LocalAA-GlobalAA-allpairs-merged-filt.tsv
zcat ${tissue}-LocalAA-GlobalAA-allpairs-merged.tsv.gz | \
  awk -F "\t" '{ if(($6 < 0.05) || ($9 < 0.05)) { print } }' \
  >> ${tissue}-LocalAA-GlobalAA-allpairs-merged-filt.tsv
gzip ${tissue}-LocalAA-GlobalAA-allpairs-merged-filt.tsv
```

## Get tied lead SNPs for each gene 

[`extract-egenes.py`](extract-egenes.py) identifies the most significant SNP(s) for each tested gene. If multiple SNPs have the same smallest p-value, all tied lead SNPs are reported. No significance threshold is applied; the lead SNP(s) is reported for **all** tested genes.  

The input should be the **unfiltered** allpairs file for a given tissue. There are two outputs:
  - `gtex.admix.global.egenes.tied.txt.gz`: Lead SNPs according to GlobalAA  
  - `gtex.admix.lava.egenes.tied.txt.gz`: Lead SNPs according to LocalAA  

## Generate some eGene sets that are repeatedly used in downstream analyses

[`get_uniq_egenes.R`](get_uniq_egenes.R) generates a few files, a few of which are specific to a particular nominal p-value cutoff:
  - `egenes_master.RData`: Lead SNPs for all tested genes; the union of `gtex.admix.global.egenes.tied.txt.gz` and `gtex.admix.lava.egenes.tied.txt.gz` for all tissues. One line per gene-SNP pair per ancestry adjustment method
  - `allp_master.RData`: Any gene-SNP pair with a nominal p-value < 0.05 (either LocalAA or GlobalAA)

Then a few files that are specific to a given nominal p-value cutoff, e.g. 1e-06:  
  - `uniq_egene_${cutoff}.RData`: For each tissue, genes that have a nominal eQTL p-value less than the cutoff for only one ancestry adjustment method. Results for all tissues are concatenated  
  - `diff_lead_${cutoff}.RData`: For each tissue, genes that have a nominal eQTL p-value less than the cutoff for BOTH ancestry adjustment methods but have different lead SNPs. Results for all tissues are concatenated  
  - `same_egene_diff_lead_${cutoff}.txt`: Unique genes from `diff_lead_${cutoff}.RData`; used to filter loci for colocalization analyses later on 

## Calculate variance in gene expression explained by local and global ancestry 
First, use [`interpolate-tss-la.py`](interpolate-tss-la.py) to extract local ancestry calls for 117AX at annotated hg19 TSSs. Then use [`variance_explained_ancestry_gene_expression.R`](variance_explained_ancestry_gene_expression.R) to calculate variance in gene expression explained by each type of ancestry, independent of the other type. This is done using a two-step regression for each gene: first, 117AX gene expression is residualized by one ancestry type (local or global). Then a second linear regression is fit on the residuals, where the predictors are the opposite ancestry type. Variance explained in the expression of a gene is taken to be the R-squared value from the second fit. This process is performed for all expressed genes (all genes in GTEx v8 normalized expression files) for each tissue. 

