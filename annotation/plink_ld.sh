#!/bin/bash

infile=$1

vcfdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno
outdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/20191030

cd ${outdir}

# # must convert VCF to bed/bim/fam for parallel computing
# for chrom in {1..22}; do  
# 	vcf=${vcfdir}/gtex.admixed.MAC10.phased.chr${chrom}.recode.vcf.gz
# 	plink --vcf ${vcf} --make-bed --out gtex.admixed.MAC10.phased.chr${chrom}
# done

# # split
# split -d -n l/15 pairs_in.txt pairs_in_

while IFS= read -r line; do
	
	var1=$(echo $line | cut -f1 -d' ')
	var2=$(echo $line | cut -f2 -d' ')

	chrom=$(echo ${var1} | cut -f1 -d'_')
	bed=plink_files/gtex.admixed.MAC10.phased.${chrom}

	d=$(plink --bfile ${bed} --ld ${var1} ${var2} --threads 2 | grep -m1 "R-sq" | sed "s/.*= //")
	if [ -z "$d" ]; then
		plink --bfile ${bed} --ld ${var1} ${var2} --threads 2 
	fi

	echo "${var1}	${var2}	${d}" >> $(basename ${infile})-out.txt
	echo "${var1} ${var2} ${d}"

done < ${infile}

