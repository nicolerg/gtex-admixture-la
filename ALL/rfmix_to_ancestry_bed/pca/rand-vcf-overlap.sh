#!/bin/bash

# For each chromosome, select a random subset of overlapping SNPs
# Used output to generate global PCs

snpdir="/mnt/lab_data/montgomery/nicolerg/rfmix-inputs"
vcfdir="/mnt/lab_data/montgomery/nicolerg/filtered"
outdir="/mnt/lab_data/montgomery/nicolerg/pca-subset"
gtex="/mnt/lab_data/montgomery/nicolerg/gtex-hg19/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_hg19.vcf.gz"

for chr in {1..22}
do 
	rand="${outdir}/rand-${chr}.txt"
	temp="${outdir}/temp-${chr}.txt"

	echo "Working on chromosome ${chr}..."

	snpin="${snpdir}/snp-bed-chr${chr}_cM.bed"
	cat ${snpin} | perl -ne 'print if (rand() < .01)' > ${rand}
	# take just the position column
	cut -f 2 ${rand} > ${temp}

	tgout="${outdir}/tg.chr${chr}.sample.vcf"
	gtexout="${outdir}/gtex.chr${chr}.sample.vcf"

	zcat ${vcfdir}/gtex.chr${chr}.vcf.gz | grep -m1 -E "^#CHROM" > ${gtexout}
	zcat ${vcfdir}/tg.chr${chr}.vcf.gz | grep -m1 -E "^#CHROM" > ${tgout}

	while read p
	do 
		tabix ${gtex} chr${chr}:${p}-${p} >> ${gtexout}
		tabix ${vcfdir}/tg.chr${chr}.vcf.gz ${chr}:${p}-${p} >> ${tgout}
		#zcat ${vcfdir}/gtex.chr$1.vcf.gz | grep -P -m1 "^chr$1\t${p}\t" >> ${gtexout}
		#zcat ${vcfdir}/tg.chr$1.vcf.gz | grep -P -m1 "^$1\t${p}\t" >> ${tgout}
	done < ${temp}

	rm ${rand}
	rm ${temp}

	bgzip ${gtexout}
	bgzip ${tgout}
done 

# now concatenate into single files
gtexmaster="${outdir}/gtex.master.subset.vcf"
tgmaster="${outdir}/tg.master.subset.vcf"

zcat ${outdir}/gtex.chr1.sample.vcf.gz > ${gtexmaster}
zcat ${outdir}/tg.chr1.sample.vcf.gz > ${tgmaster}

for chr in {2..22}
do 
	zcat ${outdir}/gtex.chr${chr}.sample.vcf.gz | tail -n +2 >> ${gtexmaster}
	zcat ${outdir}/tg.chr${chr}.sample.vcf.gz | tail -n +2 >> ${tgmaster}
done	

bgzip ${gtexmaster}

# now remove multiallelic SNPs from the 1000 Genomes intersection because tabix doesn't always return a single line for the query
sed '/|2/d' ${tgmaster} | sed '/2|/d' > ${outdir}/tg.filt.vcf
rm ${tgmaster}
mv ${outdir}/tg.filt.vcf ${outdir}/tg.master.subset.vcf
bgzip ${outdir}/tg.master.subset.vcf