# Perform local ancestry inference  

All of these steps are included in [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh).  

## Merge GTEx and 1000 Genomes VCFs

The SHAPEIT2-phased, MAF 0.01-filtered GTEx v8 WGS VCF was converted to hg19 using LiftOver for compatibility with the GRCh37 HapMap genetic map and hg19 1000 Genomes ([`liftOver-hg38Tohg19.sh`](liftOver-hg38Tohg19.sh)). Both 1000 Genomes (hg19) and GTEx v8 (hg19) VCFs were split into chromosome-level VCFs. GTEx and 1000 Genomes VCFs were then merged for each chromosome. For the manuscript, this was done with [format-rfmix.R](format-rfmix.R) (see [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh)). 

However, this merging step can be streamlined by first applying a MAF filter to the files before merging and using `bcftools` to merge the two projects. Here is an example of code that can be used to do this on a system with 22 cores. A cleaner, more efficiently parallelized version in the form of a [Snakemake pipeline](https://snakemake.readthedocs.io/en/stable/) is available [here](../admixture/Snakefile). 
```{bash}
## ALL genotyped individuals from both projects
## MAF filtered (0.05 in each project before merging)
## biallelic SNPs
## only SNPs present in both projects after MAF-filtering
## chromosome names prepended with "chr"
## hg19

outdir=/projects/nicolerg/gtex-tg-merged

# 1000 Genomes VCFs, already split by chromosome:
tg=/projects/1KG/hg19

# GTEx VCF, needs to be split by chromosome:
gtex_vcf=/projects/gtex/gtex.phased.MAF01.hg19.vcf.gz

# split GTEx VCF by chromosome (only need austosomes) and apply MAF filter
for chr in {1..22}; do
    bcftools view \
    	-r chr${chr} \
    	--min-af 0.05 \
    	-Oz \
    	-o ${outdir}/gtex_phased_chr${chr}_MAF05.vcf.gz \
    	${gtex_vcf} &
done
wait

for chr in {1..22}; do
    echo "${chr} chr${chr}" >> chr_name_conv.txt
done
format_1kg () {
	local chr=$1
	# apply MAF filter to 1000 Genomes
    bcftools view \
    	--min-af 0.05 \
    	-Oz \
    	-o ${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz \
    	${tg}/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz 
    # change 1000 Genomes naming convention
    bcftools annotate \
    	--rename-chrs chr_name_conv.txt \
    	-Oz \
    	-o ${outdir}/tg_phased_chr${chr}_MAF05_renamed.vcf.gz \
    	${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz 
}
# this will use 22 cores 
for chr in {1..22}; do format_1kg "${chr}" & done
wait

# index all VCFs
rm ${outdir}/tg_phased_chr{1..22}_MAF05.vcf.gz
for file in `ls | grep "gz"`; do
    tabix -p vcf $file 
done

# find intersection; merge result
merge_chrom () {
	mkdir -p ${outdir}/chr${chr}
    bcftools isec \
    	--collapse none \
    	-Oz \
    	-p ${outdir}/chr${chr} \
    	${outdir}/gtex_phased_chr${chr}.vcf.gz \
    	${outdir}/tg_phased_chr${chr}_MAF05_renamed.vcf.gz 
    bcftools merge \
    	-0 \
    	-m none \
    	-Oz \
    	-o ${outdir}/merged_hg38_chr${chr}.vcf.gz \
    	${outdir}/chr${chr}/0002.vcf.gz \
    	${outdir}/chr${chr}/0003.vcf.gz 
}
# this will use 22 cores 
for chr in {1..22}; do merge_chrom "${chr}" & done
```

### Generate `classes` and `alleles` RFMix inputs

[RFMix v1.5.4](https://sites.google.com/site/rfmixlocalancestryinference/) requires specially formatted inputs, including SNP positions in centimorgans. See the [Manual](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip?file_subpath=%2FRFMix_v1.5.4%2FManual.pdf) for more details. While the merging steps of [format-rfmix.R](format-rfmix.R) can be improved (see above), it is still helpful for generating these inputs.   

[format-rfmix.R](format-rfmix.R) requires several files:  
  - [gtex-admixed0.9.txt](../metadata/gtex-admixed0.9.txt): list of 117AX GTEx IDs
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

Outputs of [format-rfmix.R](format-rfmix.R):
  - `alleles-chr${CHROM}`: One row per SNP and one column per haplotype (see `alleles` input spec in the [RFMix Manual](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip?file_subpath=%2FRFMix_v1.5.4%2FManual.pdf))
  - `snp-bed-chr${CHROM}` where the columns are `chromosome`, `start (bp)`, `stop (bp)`, `rsid` (hg19 coordinates), e.g.:
    ```
    chr10 60753 60754 rs554788161
    chr10 64774 64775 rs368644811
    chr10 64851 64852 rs548882871
    ```
  - `classes-chr${CHROM}`: One row with one column per haplotype (see `classes` input spec in the [RFMix Manual](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip?file_subpath=%2FRFMix_v1.5.4%2FManual.pdf))
  - `subject-ids-chr${CHROM}`: One row per GTEx ID, in the order of the columns in `alleles-chr${CHROM}`

`classes-chr${CHROM}` and `subject-ids-chr${CHROM}` are the same for all chromosomes, so these outputs can be condensed into a single file once all chromosomes are processed, as in [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh):
```
# condense class and subject ID files
mv ${outdir}/$(ls ${outdir} | grep "classes-chr" | head -1) ${outdir}/classes
rm -f ${outdir}/classes-chr*
mv ${outdir}/$(ls ${outdir} | grep "subject-ids-chr" | head -1) ${outdir}/subject-ids
# mv ${outdir}/subject-ids-chr1 ${outdir}/subject-ids
rm -f ${outdir}/subject-ids-chr*
```

## Convert SNP positions from bp to cM to generate `snp_locations` RFMix input 

Next, SNP positions need to be converted from bp to centimorgans (cM). First, download a genetic map. This manuscript used the GRCh37 HapMap genetic map, available [here](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/).  

**Note:** A GRCh38 genetic map can now be downloaded from [Beagle](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).  

Follow the code block in [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh) labelled **convert bp to genetic distance**, which calls [bp-to-cM.py](bp-to-cM.py).  

**Note:** [bp-to-cM.py](bp-to-cM.py) is a modified version of code written by Joe Pickrell, available [here](https://github.com/joepickrell/1000-genomes-genetic-maps/blob/master/scripts/interpolate_maps.py).

## Run RFMix

Download [RFMix v1.5.4](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip). Then follow the code block in [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh) labelled **run RFMix**. 

**Note:** When RFMix was originally run for this manuscript, many more SNPs were included than necessary. This required RFMix to be parallelized per individual. The code in [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh) reflects this. Alternatively, apply a MAF filter and LD-prune your VCFs to get \~1M total SNPs before running RFMix, in which case you do not need to split up the analysis for each individual. It's fine to include more SNPs (RFMix is not affected by LD since SNP positions are modelled in cM), but it certainly does slow it down.   

## Collapse local ancestry calls into BED files of contiguous local ancestry 

Generate 2 BED files per individuals (1 per haplotype). Follow the code chunk in [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh) that calls [collapse_ancestry_single_chr.py](collapse_ancestry_single_chr.py).  

**Note:** [collapse_ancestry_single_chr.py](collapse_ancestry_single_chr.py) is a modified version of code written by Alicia R. Martin, available [here](https://github.com/armartin/ancestry_pipeline/blob/master/collapse_ancestry.py).

## Calculate global ancestry proportions  

Scripts from Alicia R. Martin were also used to [calculate global ancestry](https://github.com/armartin/ancestry_pipeline/blob/master/lai_global.py) and [plot karyograms](https://github.com/armartin/ancestry_pipeline/blob/master/plot_karyogram.py) for each individual. See her repository for other scripts related to processing RFMix outputs: https://github.com/armartin/ancestry_pipeline.  

## Make master local ancestry files

See the last chunk of [format-rfmix-pipeline.sh](format-rfmix-pipeline.sh) to generate chromosome-level local ancestry BED files.   

Call [interpolate-local-anc.py](interpolate-local-anc.py) to extract local ancestry for each SNP in a target VCF (`subset` path). `mapfile` is each chromosome-level local ancestry BED file. See how target VCFs were generated for this manuscript [here](#prepare-vcfs). 
