#!/bin/bash

set -e 

# get 1000 Genomes European, Eastern Asian, and African allele frequencies for all SNPs tested in the LAVA pipelines

outdir=$1

all_snps=$outdir/all_rsid.txt 
zcat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.rsid.map.txt.gz | sed -e '1d' | cut -f1 > $outdir/all_rsid.txt 

kg=/mnt/lab_data/montgomery/shared/1KG/hg38

for chrom in {1..22}; do

	vcf=${kg}/ALL.chr${chrom}_GRCh38.genotypes.20170504.vcf.gz
	outfile=${outdir}/gtex.admixed.chr${chrom}.1kg.AF

	vcftools --gzvcf ${vcf} --out ${outfile} --snps ${all_snps} --get-INFO EAS_AF --get-INFO EUR_AF --get-INFO AFR_AF &

done

wait

echo "CHROM	POS	REF	ALT	EAS_AF	EUR_AF	AFR_AF" > ${outdir}/gtex.admixed.1kg.hg38.AF.INFO
for chrom in {1..22}; do 

	sed -e '1d' ${outdir}/gtex.admixed.chr${chrom}.1kg.AF.INFO >> ${outdir}/gtex.admixed.1kg.hg38.AF.INFO

done

# now match to GTEx variant_id
python ~/gtex-admix/scripts/LAVA/annotation/match_af_gtex_id.py $outdir 
