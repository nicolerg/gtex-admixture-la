#!/bin/bash

cd /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/admixture_pruned/gtex_hg38
infile=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz

for chr in {1..22}; do 
	zgrep -E -m1 "^#CHROM" ${infile} > gtex.hg38.phased.chr${chr}.vcf
done

for chr in {1..22}; do 
	zgrep -E "^chr${chr}	" ${infile} >> gtex.hg38.phased.chr${chr}.vcf &
done
wait 
