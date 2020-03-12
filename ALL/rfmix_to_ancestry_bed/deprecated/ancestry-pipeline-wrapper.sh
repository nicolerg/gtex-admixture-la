#!/bin/bash

# Author: Nicole Gay 
# Winter 2018
# Process outputs from RFMix using scripts from https://github.com/armartin/ancestry_pipeline
# 1. Generate 2 bed files per individual (1 per haplotype) with ancestry assignments
## 2. Plot karyotypes for all individuals (plot_karyogram.py from armartin, with edits from me)
# 3. Calculate global ancestry for all indivuals 

base="/mnt/lab_data/montgomery/nicolerg/rfmix-admixed"
outdir="${base}/bed"
indir="${base}/whole-outputs"
ids="${base}/gtex-ids"
srcdir="/users/nicolerg/gtex-admix/scripts/ancestry_pipeline"

# while read ind
# do
# 	echo "Working on ${ind}..."
# 	echo "Collapsing results into .bed files"
# 	# generate _A.bed and _B.bed for every individual
# 	taskset -c 24-47 python ${srcdir}/collapse_ancestry.py \
# 	--rfmix "${indir}/admixed-chr1.0.Viterbi.txt.gz" \
# 	--snp_locations "${base}/inputs/snp_locations-chr1.gz" \
# 	--ind ${ind} \
# 	--ind_info "${base}/inputs/subject-ids" \
# 	--pop_labels 'ASN,EUR,AFR' \
# 	--out "${outdir}/${ind}" \
# 	--fbk "${indir}/admixed-chr1.0.ForwardBackward.txt.gz" \
# 	--fbk_threshold 0.9 &
# done < ${ids}

# wait

# # while read ind
# # do
# # 	echo "Plotting karyotypes"
# # 	# plot karyotypes
# # 	python ${srcdir}/plot_karyogram.py \
# # 	--bed_a "${outdir}/${ind}_A.bed" \
# # 	--bed_b "${outdir}/${ind}_B.bed" \
# # 	--ind ${ind} \
# # 	--centromeres "${srcdir}/centromeres_hg19.bed" \
# # 	--pop_order 'ASN,EUR,AFR' \
# # 	--out "${outdir}/${ind}-karyo.png" &
# # 	echo ""
# # done < ${ids}

# # manually make list of bed files:
# # two file names per line, corresponding to haps A and B for the same ind --> gtex.bed.list.txt

# echo "Calculating global ancestry"
# python ${srcdir}/lai_global.py \
# --bed_list ${base}/bed/gtex.bed.list.txt \
# --ind_list ${ids} \
# --pops 'ASN,EUR,AFR,UNK' \
# --out ${outdir}/gtex_ai_global.txt

# # plot haplotypes
# srcdir="/users/nicolerg/gtex-admix/scripts/rfmix_analysis"
# Rscript ${srcdir}/plot-haplotypes.R

# # combine bed files into master ancestry files
# for file in `ls ${outdir} | grep "GTEX-"`; do
# 	ind=`echo ${file} | sed "s/\..*//"`
# 	cut -f 1,2,3,4 ${outdir}/${file} | sed "s/$/	${ind}/" | sed "s/^/chr/" > ${outdir}/master/${file}.tmp
# done

# for chr in {1..22}; do
# 	outfile=${outdir}/master/ancestry.chr${chr}.bed
# 	for file in `ls ${outdir}/master | grep "tmp"`; do
# 		grep -E "^chr${chr}	" ${outdir}/master/${file} >> ${outfile}
# 	done
# done

# rm ${outdir}/master/*bed.tmp

# cp /mnt/lab_data/montgomery/nicolerg/rfmix-outputs/haplotype-bed/ancestry.chr* /mnt/lab_data/montgomery/nicolerg/admixed-ancestry

# for file in `ls ${outdir}/master | grep "ancestry.chr"`; do 
# 	cat ${outdir}/master/${file} >> /mnt/lab_data/montgomery/nicolerg/admixed-ancestry/${file}
# done

# wait

# make master ancestry files from haplotype bed files
for chr in {1..22}; do
	taskset -c 24-47 python /users/nicolerg/gtex-admix/scripts/LAVA/interpolate-local-anc.py ${chr} &

	# running=`ps -ef | grep "nicolerg" | grep "interpolate-local-anc" | wc -l`
	# while [ $running -gt 5 ]; do
	# 	sleep 30
	# 	running=`ps -ef | grep "nicolerg" | grep "interpolate-local-anc" | wc -l`
	# done
done
