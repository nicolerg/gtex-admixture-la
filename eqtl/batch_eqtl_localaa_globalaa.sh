#!/bin/bash

# run eQTL calling for a single tissue with both LocalAA and GlobalAA 

set -e 

base=~/local-eqtl/admixed # root data folder
srcdir=~/gtex-admix/scripts/LAVA
vcfdir=${base}/geno # filtered down to 117AX; MAC10 filter
exprdir=GTEx_Analysis_v8_eQTL_expression_matrices # GTEx v8 normalized expression for each tissue
admix_ids=gtex-admixed0.9.txt # see metadata folder in this repository 
covdir= # same as ${outdir} for concat_cov.R
la_cov= # directory with LA covariates

tissue=Adipose_Subcutaneous # set a tissue. must correspond to GTEx file prefix

#=================================================================

# first, get everything set up to run eQTL calling 
	
exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz
globalcov=${covdir}/${tissue}.all_covariates.txt

# set up directories
indir=${base}/${tissue}
if [ ! -d "$indir" ]; then
	mkdir $indir
fi
edir=${indir}/chunks
if [ ! -d "$edir" ]; then
	mkdir $edir
fi
for chr in {1..22}; do
	outdir=${indir}/chr${chr}
	if [ ! -d "$outdir" ]; then
		mkdir $outdir
	fi
done

# split up expression
for chr in {1..22}; do
	geno=${vcfdir}/gtex.admixed.MAC10.phased.chr${chr}.recode.vcf.gz
	outdir=${indir}/chr${chr}

	# split expression into chunks
	if [ ! -e "${indir}/chunks/${tissue}.chr${chr}.expression.03.gz" ]; then 
		exp=${indir}/${tissue}.chr${chr}.expression.bed
		zgrep -E "^chr${chr}	" ${exprfile} > ${exp}
		split -d -l 50 ${exp} "${edir}/${tissue}.chr${chr}.expression."

		header=${edir}/header.txt
		zcat ${exprfile} | head -1 > ${header}
		for file in `ls ${edir} | grep "chr${chr}.expression"`
		do
			cat ${header} > ${edir}/tmp.chr${chr}.txt
			cat ${edir}/${file} >> ${edir}/tmp.chr${chr}.txt
			rm ${edir}/${file}
			mv ${edir}/tmp.chr${chr}.txt ${edir}/${file}
		done
		gzip ${edir}/${tissue}.chr${chr}.expression*
	fi
done

#=================================================================

# now, run eQTL calling 
	
for chr in {1..22}; do

	outdir=${indir}/chr${chr}
	geno=${vcfdir}/gtex.admixed.MAC10.phased.chr${chr}.recode.vcf.gz
	localcov=${la_cov}/chr${chr}.hg19.maf0.localcov.tsv.gz

	for chunk in `ls ${edir} | grep "chr${chr}.expression"`; do
		gene=${edir}/${chunk}
		suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`
		echo "Starting chr $chr chunk $suf"
		out="${outdir}/${tissue}-LocalAA-GlobalAA-allpairs-chr${chr}-${suf}.tsv"
		taskset -c 38-47 Rscript ${srcdir}/eqtl_localaa_globalaa.R ${chr} ${tissue} ${gene} ${globalcov} ${geno} ${localcov} ${out} &

		# this is just to limit the number of jobs launched on a cluster without a job queue 
		running=`ps -ef | grep "$USER" | grep "eqtl_localaa_globalaa.R" | wc -l`
		while [ $running -gt 9 ]; do
			sleep 180
			running=`ps -ef | grep "$USER" | grep "eqtl_localaa_globalaa.R" | wc -l`
		done
		
	done
done

wait

#=================================================================

# merge chunks (per chromosome)

for chr in {1..22}; do
	outdir=${indir}/chr${chr}
	f=`ls ${outdir} | grep "${tissue}-LocalAA-GlobalAA-allpairs" | head -1`
	head -1 ${outdir}/${f} > ${outdir}/${tissue}-LocalAA-GlobalAA-merged-allpairs-chr${chr}.tsv
	for chunk in `ls ${outdir} | grep "${tissue}-LocalAA-GlobalAA-allpairs-chr${chr}"`; do
		file=${outdir}/${chunk}
		cat ${file} | sed -e '1d' >> ${outdir}/${tissue}-LocalAA-GlobalAA-merged-allpairs-chr${chr}.tsv
		#rm ${file}
	done
	gzip ${outdir}/${tissue}-LocalAA-GlobalAA-merged-allpairs-chr${chr}.tsv
done

# #=================================================================

# merge into whole genome

outdir=${indir}

zcat ${outdir}/chr1/${tissue}-LocalAA-GlobalAA-merged-allpairs-chr1.tsv.gz | head -1 > ${outdir}/${tissue}-LocalAA-GlobalAA-merged-allpairs.tsv
for dir in chr{1..22}; do
	file=${outdir}/${dir}/${tissue}-LocalAA-GlobalAA-merged-allpairs-${dir}.tsv.gz
	zcat $file | sed -e '1d' >> ${outdir}/${tissue}-LocalAA-GlobalAA-merged-allpairs.tsv
done
gzip ${outdir}/${tissue}-LocalAA-GlobalAA-merged-allpairs.tsv

