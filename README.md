# Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx

This repository provides code and resources used to perform the analyses presented in the bioRxiv preprint [Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx](https://www.biorxiv.org/content/10.1101/836825v1) by Nicole Gay, et al. The author of all code is Nicole Gay (nicolerg@stanford.edu) unless otherwise noted. Please see the manuscript for context.   

Note that code is primarily provided to replicate analyses in the manuscript. However, when applicable, updated or more efficient versions of the code are also provided.  

"117AX" refers to a subset of 117 admixed invididuals in the GTEx v8 release whose GTEx IDs are provided in [`metadata/gtex-admix0.9.txt`](metadata/gtex-admix0.9.txt). These scripts could be adapted for any subset of individuals.  

## Table of contents  
- **Perform local ancestry inference:** Use RFMix v1.5.4 and 1000 Genomes reference panels to perform local ancestry inference on 117AX   
- **eQTL calling:** Perform eQTL calling with both Global Ancestry Adjustment (GlobalAA) and Local Ancestry Adjustment (LocalAA)  
- **Colocalization:** <Add something here>

## Perform local ancestry inference  

All of these steps are included in [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh).  

### Merge GTEx and 1000 Genomes VCFs

The SHAPEIT2-phased, MAF 0.01-filtered GTEx v8 WGS VCF was converted to hg19 using LiftOver for compatibility with the GRCh37 HapMap genetic map and hg19 1000 Genomes ([`liftOver-hg38Tohg19.sh`](rfmix/liftOver-hg38Tohg19.sh)). Both 1000 Genomes (hg19) and GTEx v8 (hg19) VCFs were split into chromosome-level VCFs. GTEx and 1000 Genomes VCFs were then merged for each chromosome. For the manuscript, this was done with [format-rfmix.R](rfmix/format-rfmix.R) (see [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh)). 

However, this merging step can be streamlined by first applying a MAF filter to the files before merging and using `bcftools` to merge the two projects. Here is an example of code that can be used to do this on a system with 22 cores. A cleaner, more efficiently parallelized version in the form of a [Snakemake pipeline](https://snakemake.readthedocs.io/en/stable/) is available [here](rfmix/Snakefile). 
```
## ALL genotyped individuals from both projects
## MAF filtered (0.05 in each project before merging)
## biallelic SNPs
## only SNPs present in both projects after MAF-filtering
## chromosome names prepended with "chr"
## hg19

outdir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged

# 1000 Genomes VCFs, already split by chromosome:
tg=/mnt/lab_data/montgomery/shared/1KG/hg19

# GTEx VCF, needs to be split by chromosome:
gtex_vcf=/mnt/lab_data/montgomery/nicolerg/gtex.phased.MAF01.hg19.vcf.gz

# split GTEx VCF by chromosome (only need austosomes) and apply MAF filter
for chr in {1..22}; do
    bcftools view -r chr${chr} --min-af 0.05 -Oz -o ${outdir}/gtex_phased_chr${chr}_MAF05.vcf.gz ${gtex_vcf} &
done
wait

# apply MAF filter to 1000 Genomes
for chr in {1..22}; do
    bcftools view --min-af 0.05 -Oz -o ${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz ${tg}/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz &
done
wait

# change 1000 Genomes naming convention
for chr in {1..22}; do
    echo "${chr} chr${chr}" >> chr_name_conv.txt
done
for chr in {1..22}; do
    bcftools annotate --rename-chrs chr_name_conv.txt -Oz -o ${outdir}/tg_phased_chr${chr}_MAF05_renamed.vcf.gz ${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz &
done
wait

# index all VCFs
rm ${outdir}/tg_phased_chr{1..22}_MAF05.vcf.gz
for file in `ls | grep "gz"`; do
    tabix -p vcf $file &
done
wait

# find intersection; merge result
for chr in {1..22}; do
    mkdir -p ${outdir}/chr${chr}
done
for chr in {1..22}; do
    bcftools isec --collapse none -Oz -p ${outdir}/chr${chr} ${outdir}/gtex_phased_chr${chr}.vcf.gz ${outdir}/tg_phased_chr${chr}_MAF05_renamed.vcf.gz &
done
wait
for chr in {1..22}; do
    bcftools merge -0 -m none -Oz -o ${outdir}/merged_hg38_chr${chr}.vcf.gz ${outdir}/chr${chr}/0002.vcf.gz ${outdir}/chr${chr}/0003.vcf.gz &
done
```

### Generate `classes` and `alleles` RFMix inputs

[RFMix v1.5.4](https://sites.google.com/site/rfmixlocalancestryinference/) requires specially formatted inputs, including SNP positions in centimorgans. See the [Manual](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip?file_subpath=%2FRFMix_v1.5.4%2FManual.pdf) for more details. While the merging steps of [format-rfmix.R](rfmix/format-rfmix.R) can be improved (see above), it is still helpful for generating these inputs.   

[format-rfmix.R](rfmix/format-rfmix.R) requires several files:  
  - (gtex-admixed0.9.txt)[metadata/gtex-admixed0.9.txt]: list of 117AX GTEx IDs
  - 1000 Genomes metadata file with population labels, i.e.:
    ```
    Subject Subpop  Superpop
    HG02308 ACB AFR
    HG02255 ACB AFR
    ...
    ```
  - `CHROM`, number to specify the chromosome 
  - GTEx VCF (hg19) for a single chromsome (`CHROM`)
  - 1000 Genomes VCF (hg19) for a single chromosome (`CHROM`)

Outputs of [format-rfmix.R](rfmix/format-rfmix.R):
  - `alleles-chr${CHROM}`: One row per SNP and one column per haplotype (see `alleles` input spec in the [RFMix Manual](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip?file_subpath=%2FRFMix_v1.5.4%2FManual.pdf))
  - `snp-bed-chr${CHROM}` where the columns are `chromosome`, `start (bp)`, `stop (bp)`, `rsid` (hg19 coordinates), e.g.:
    ```
    chr10 60753 60754 rs554788161
    chr10 64774 64775 rs368644811
    chr10 64851 64852 rs548882871
    ```
  - `classes-chr${CHROM}`: One row with one column per haplotype (see `classes` input spec in the [RFMix Manual](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip?file_subpath=%2FRFMix_v1.5.4%2FManual.pdf))
  - `subject-ids-chr${CHROM}`: One row per GTEx ID, in the order of the columns in `alleles-chr${CHROM}`

`classes-chr${CHROM}` and `subject-ids-chr${CHROM}` are the same for all chromosomes, so these outputs can be condensed into a single file once all chromosomes are processed, as in [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh):
```
# condense class and subject ID files
mv ${outdir}/$(ls ${outdir} | grep "classes-chr" | head -1) ${outdir}/classes
rm -f ${outdir}/classes-chr*
mv ${outdir}/$(ls ${outdir} | grep "subject-ids-chr" | head -1) ${outdir}/subject-ids
# mv ${outdir}/subject-ids-chr1 ${outdir}/subject-ids
rm -f ${outdir}/subject-ids-chr*
```

### Convert SNP positions from bp to cM to generate `snp_locations` RFMix input 

Next, SNP positions need to be converted from bp to centimorgans (cM). First, download a genetic map. This manuscript used the GRCh37 HapMap genetic map, available [here](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/).  

**Note:** A GRCh38 genetic map can now be downloaded from [Beagle](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).  

Follow the code block in [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh) labelled **convert bp to genetic distance**, which calls [bp-to-cM.py](rfmix/bp-to-cM.py).  

**Note:** [bp-to-cM.py](rfmix/bp-to-cM.py) is a modified version of code written by Joe Pickrell, available [here](https://github.com/joepickrell/1000-genomes-genetic-maps/blob/master/scripts/interpolate_maps.py).

### Run RFMix

Download [RFMix v1.5.4](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip). Then follow the code block in [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh) labelled **run RFMix**. 

**Note:** When RFMix was originally run for this manuscript, many more SNPs were included than necessary. This required RFMix to be parallelized per individual. The code in [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh) reflects this. Alternatively, apply a MAF filter and LD-prune your VCFs to get \~1M total SNPs before running RFMix, in which case you do not need to split up the analysis for each individual. It's fine to include more SNPs (RFMix is not affected by LD since SNP positions are modelled in cM), but it certainly does slow it down.   

### Collapse local ancestry calls into BED files of contiguous local ancestry 

Generate 2 BED files per individuals (1 per haplotype). Follow the code chunk in [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh) that calls [collapse_ancestry_single_chr.py](rfmix/collapse_ancestry_single_chr.py).  

**Note:** [collapse_ancestry_single_chr.py](rfmix/collapse_ancestry_single_chr.py) is a modified version of code written by Alicia R. Martin, available [here](https://github.com/armartin/ancestry_pipeline/blob/master/collapse_ancestry.py).

### Calculate global ancestry proportions  

Scripts from Alicia R. Martin were also used to [calculate global ancestry](https://github.com/armartin/ancestry_pipeline/blob/master/lai_global.py) and [plot karyograms](https://github.com/armartin/ancestry_pipeline/blob/master/plot_karyogram.py) for each individual. See her repository for other scripts related to processing RFMix outputs: https://github.com/armartin/ancestry_pipeline.  

### Make master local ancestry files

See the last chunk of [format-rfmix-pipeline.sh](rfmix/format-rfmix-pipeline.sh) to generate chromosome-level local ancestry BED files.   

Call [interpolate-local-anc.py](rfmix/interpolate-local-anc.py) to extract local ancestry for each SNP in a target VCF (`subset` path). `mapfile` is each chromosome-level local ancestry BED file. See how target VCFs were generated for this manuscript [here](#prepare-vcfs). 

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

### Get tied lead SNPs for each gene 

[`extract-egenes.py`](eqtl/extract-egenes.py) identifies the most significant SNP(s) for each tested gene. If multiple SNPs have the same smallest p-value, all tied lead SNPs are reported. No significance threshold is applied; the lead SNP(s) is reported for **all** tested genes.  

There are two outputs:
  - `gtex.admix.global.egenes.tied.txt.gz`: Lead SNPs according to GlobalAA  
  - `gtex.admix.lava.egenes.tied.txt.gz`: Lead SNPs according to LocalAA  





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