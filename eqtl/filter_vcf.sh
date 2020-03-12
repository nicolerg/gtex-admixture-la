#!/bin/bash

# filter VCFs by AC (only do once per group of subjects)

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

