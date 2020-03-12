#!/bin/bash
# Outputs one headerless VCF and corresponding index file for each chromosome for both GTEx (hg19) and 1000 Genomes (b37)

outdir=/mnt/lab_data/montgomery/nicolerg/filtered
gtex=/mnt/lab_data/montgomery/nicolerg/gtex-hg19/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_hg19.vcf.gz
tgbase=/mnt/lab_data/montgomery/shared/1KG

for chr in {1..22}
do 
	echo "Working on chr ${chr}"
	tout=${outdir}/tg.chr${chr}.vcf
	gout=${outdir}/gtex.chr${chr}.vcf
	tgin=${tgbase}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	zcat ${gtex} | grep -P "^#CHROM|^chr${chr}\t" > ${gout}
	bgzip ${gout}
	tabix -fp vcf ${outdir}/gtex.chr${chr}.vcf.gz
	zcat ${tgin} | grep -v '##' > ${tout}
	bgzip ${tout}
	tabix -fp vcf ${outdir}/tg.chr${chr}.vcf.gz
done 
