#!/bin/bash
# LiftoverVcf for all 1000 Genomes chromosomes

b37="/mnt/lab_data/montgomery/shared/1KG"
refhg38="/mnt/lab_data/montgomery/shared/genomes/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
chain19to38="/mnt/lab_data/montgomery/shared/liftOver/chains/hg19ToHg38.over.chain.gz"
indir="/mnt/lab_data/montgomery/nicolerg/1KG-hg19"
zcat ${chain19to38} > "${indir}/hg19ToHg38.over.chain"
outdir="/mnt/lab_data/montgomery/nicolerg/1KG-hg38"

# for all chromosomes
for n in {1..22}
do

	if [ $n -eq 6 ]
	then

		input="/mnt/lab_data/montgomery/nicolerg/1KG-hg19/chr6.b37.vcf.gz"

	else

		input="${b37}/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

	fi

	hg19="${indir}/chr${n}.formatted.vcf"
	zcat ${input} | sed "s/^${n}/chr${n}/" > "${hg19}"
	output="${outdir}/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a_hg38.20130502.genotypes.vcf"

	java -jar /users/nicolerg/picard.jar LiftoverVcf \
	     I="${hg19}" \
	     O="${output}" \
	     CHAIN="${chain19to38}" \
	     REJECT="${outdir}/chr${n}_rejected_variants.vcf" \
	     R="${refhg38}" # target reference build (hg38)

	rm ${hg19}
	gzip ${output}

done


# # for a single chromosome:
# n=21
# input="${b37}/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
# hg19="${indir}/chr${n}.formatted.vcf"
# zcat ${input} | sed "s/^${n}/chr${n}/" > "${hg19}"
# output="${outdir}/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a_hg38.20130502.genotypes.vcf"

# java -jar picard.jar LiftoverVcf \
#      I="${hg19}" \
#      O="${output}" \
#      CHAIN="${chain19to38}" \
#      REJECT="${outdir}/chr${n}_rejected_variants.vcf" \
#      R="${refhg38}" # target reference build (hg38)

# #rm ${hg19}
# gzip ${output}


# indir="/mnt/lab_data/montgomery/nicolerg/1KG-hg19"
# outdir="/mnt/lab_data/montgomery/nicolerg/1KG-hg38"
# build38="/mnt/lab_data/montgomery/shared/genomes/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
# #build38="/mnt/lab_data/montgomery/nicolerg/hg38-ref"
# #build19="/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa"
# chain19to38="/mnt/lab_data/montgomery/nicolerg/hg38-ref/hg19ToHg38.over.chain"
# #chainb37tohg19="/users/nicolerg/gtex-admix/data/b37tohg19.chain"

# # lift over directly from b37 (prepend chr) but use whole reference 

# NUMBERS="2 8 11 12 15 18 19 22"

# for n in `echo $NUMBERS`
# do 
	
# 	infile="${indir}/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
# 	input="${indir}/chr${n}.hg19.vcf"
# 	zcat ${infile} | sed "s/^${n}/chr${n}/" > "${input}"
# 	output="${outdir}/chr${n}.hg38.vcf"
# 	reference="${build38}"

# 	java -jar picard.jar LiftoverVcf \
# 	     I="${input}" \
# 	     O="${output}" \
# 	     CHAIN="${chain19to38}" \
# 	     REJECT="${outdir}/chr${n}_rejected_variants.vcf" \
# 	     R="${reference}" # target reference build (hg38)

# 	gzip ${output}

# done

# for n in {1..22}
# do 

# 	input="${indir}/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# 	output="${indir}/chr${n}.hg19.vcf"
# 	reference="${build19}"

# 	java -jar picard.jar LiftoverVcf \
# 	     I="${input}" \
# 	     O="${output}" \
# 	     CHAIN="${chainb37tohg19}" \
# 	     REJECT="${outdir}/chr${n}_rejected_variants.vcf" \
# 	     R="${reference}" # target reference build (hg38)

# 	gzip ${output}

# 	input="${indir}/chr${n}.hg19.vcf.gz"
# 	output="${outdir}/chr${n}.hg38.vcf"
# 	reference="${build38}/chr${n}.fa"

# 	java -jar picard.jar LiftoverVcf \
# 	     I="${input}" \
# 	     O="${output}" \
# 	     CHAIN="${chain19to38}" \
# 	     REJECT="${outdir}/chr${n}_pass2_rejected_variants.vcf" \
# 	     R="${reference}" # target reference build (hg38)

# 	gzip ${output}

# done

# input="${indir}/chr6.b37.vcf.gz"
# n=6

# output="${indir}/chr${n}.hg19.vcf"
# reference="${build19}"

# java -jar picard.jar LiftoverVcf \
#      I="${input}" \
#      O="${output}" \
#      CHAIN="${chainb37tohg19}" \
#      REJECT="${outdir}/chr${n}_rejected_variants.vcf" \
#      R="${reference}" # target reference build (hg38)

# gzip ${output}

# input="${indir}/chr${n}.hg19.vcf.gz"
# output="${outdir}/chr${n}.hg38.vcf"
# reference="${build38}"

# java -jar picard.jar LiftoverVcf \
#      I="${input}" \
#      O="${output}" \
#      CHAIN="${chain19to38}" \
#      REJECT="${outdir}/chr${n}_pass2_rejected_variants.vcf" \
#      R="${reference}" # target reference build (hg38)

# gzip ${output}


# NUMBERS="2 8 11 12 15 18 19 22"

# for n in `echo $NUMBERS`  
# do

# 	input="${indir}/chr${n}.hg19.vcf.gz"
# 	output="${outdir}/chr${n}.hg38.vcf"
# 	reference="${build38}"

# 	java -jar picard.jar LiftoverVcf \
# 	     I="${input}" \
# 	     O="${output}" \
# 	     CHAIN="${chain19to38}" \
# 	     REJECT="${outdir}/chr${n}_pass2_rejected_variants.vcf" \
# 	     R="${reference}" # target reference build (hg38)

# 	gzip ${output}

# done

# input="${indir}/chr6.hg19.vcf.gz"
# output="${outdir}/chr6.hg38.vcf"
# reference="${build38}"

# java -jar picard.jar LiftoverVcf \
#      I="${input}" \
#      O="${output}" \
#      CHAIN="${chain19to38}" \
#      REJECT="${outdir}/chr6_pass2_rejected_variants.vcf" \
#      R="${reference}" # target reference build (hg38)

# gzip ${output}

