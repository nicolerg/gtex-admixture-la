# Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx

This repository provides code and resources used to perform the analyses presented in the bioRxiv preprint [Impact of admixture and ancestry on eQTL analysis and GWAS colocalization in GTEx](https://www.biorxiv.org/content/10.1101/836825v1) by Nicole Gay, et al. The author of all code is Nicole Gay (nicolerg@stanford.edu) unless otherwise noted.  

Note that code is primarily provided to replicate analyses in the manuscript. However, when applicable, updated or more efficient versions of the code are also provided.  

## Table of contents  
- **Perform local ancestry inference:** Use RFMix v1.5.4 and 1000 Genomes reference panels to perform local ancestry inference on 117 admixed individuals in GTEx   
- **eQTL calling:** Perform eQTL calling with both Global Ancestry Adjustment (GlobalAA) and Local Ancestry Adjustment (LocalAA)  
- **Colocalization:** <Add something here>

## Perform local ancestry inference  

### Merge GTEx and 1000 Genomes VCFs

The SHAPEIT2-phased, MAF 0.01-filtered GTEx v8 WGS VCF was converted to hg19 using LiftOver for compatibility with the GRCh37 HapMap genetic map and hg19 1000 Genomes ([`liftOver-hg38Tohg19.sh`](rfmix/liftOver-hg38Tohg19.sh)). Both 1000 Genomes (hg19) and GTEx v8 (hg19) VCFs were split into chromosome-level VCFs. GTEx and 1000 Genomes VCFs were then merged for each chromosome. For the manuscript, this was done with [format-rfmix.R](rfmix/format-rfmix.R) (see [format-rfmix.sh](rfmix/format-rfmix.sh)). 

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

### Format inputs for RFMix

[RFMix v1.5.4](https://sites.google.com/site/rfmixlocalancestryinference/) requires specially formatted inputs, including SNP positions in centimorgans. See the [Manual](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip?file_subpath=%2FRFMix_v1.5.4%2FManual.pdf) for more details. While the merging steps of [format-rfmix.R](rfmix/format-rfmix.R) can be improved (see above), it is still helpful for generating these inputs.  

[format-rfmix.R](rfmix/format-rfmix.R) requires several files:  
  - 





The GRCh37 HapMap genetic map was downloaded from [here](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/).  

**Note:** A GRCh38 genetic map can now be downloaded from [Beagle](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

### Download RFMix v1.5.4

Download the zip file found [here](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip).

### Format inputs for RFMix
Run [`format-rfmix.sh`](rfmix/format-rfmix.sh) to do the following:
  - Format inputs for RFMix (uses [`format-rfmix.R`](rfmix/format-rfmix.R))
  - Run RFMix (v1.5.4)
  - Generate per-individual BED files of contiguous ancestry assignments (`collapse_ancestry_single_chr.py`)
  - Collapse all per-individual BED files into a master BED file (one per chromosome)
  - Use the master BED file to interpolate local ancestry for every variant specified in the VCF file for your admixed individuals (per chromosome). This file is used to construct the local ancestry variables in eQTL calling (`interpolate-local-anc.py`)


## eQTL calling (eqtl subdirectory)

### Prepare VCFs 

The SHAPEIT2-phased, MAF 0.01-filtered GTEx v8 WGS VCF was converted to hg19 using LiftOver for compatibility with the hg19 HapMap genetic map [`liftOver-hg38Tohg19.sh`](rfmix/liftOver-hg38Tohg19.sh). This VCF was then split into chromosome-level VCFs, each of which was then filtered to include only 117AX samples and MAC > 10 using `vcftools`:
```
vcfin=genotypes-hg19
vcfdir=admixed/geno

for chr in {1..22}; do
	vcf_unfilt=${vcfin}/gtex.MAF01.phased.chr${chr}.vcf.gz
	vcf_out=${vcfdir}/gtex.admixed.MAC10.phased.chr${chr}
	taskset -c 24-47 vcftools \
		--gzvcf ${vcf_unfilt} \
		--keep ${admix_ids} \
		--recode-INFO-all \
		--mac 10 \
		--max-alleles 2 \
		--min-alleles 2 \
		--out ${vcf_out} \
		--recode &
done

wait

gzip ${vcfdir}/gtex.admixed.MAC10.phased.chr*
```

### Prepare covariates

[`prepare_covariates.sh`](eqtl/prepare_covariates.sh) concatenates covariates provided by GTEx and new PEER factors. It requires paths to two GTEx v8 files:  
  - `${tissue}.v8.normalized_expression.bed.gz`
  - `${tissue}.v8.covariates.txt`  
It also calls 3 other scripts:  
  - [`filter_expression_admixed.R`](eqtl/filter_expression_admixed.R)
  - [`concat_cov.R`](eqtl/concat_cov.R)
  - [`run_PEER.R`](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R) 

### Run eQTL calling with both LocalAA and GlobalAA 

`batch_eqtl_localaa_globalaa.sh`
calls:
  - `eqtl_localaa_globalaa.R`

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