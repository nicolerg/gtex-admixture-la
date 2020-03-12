#!/bin/bash

set -e 

################################################################################################################
## Run RFMix
################################################################################################################

outdir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_outputs
indir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_inputs
subj=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/gtex-admixed0.9.txt

mkdir -p ${outdir}/whole-outputs

# # this block runs the algorithm separately for each admixed individual
# # it iteratively parses the input files to include only one admixed individual and all reference panels

# # just do this for the first 10 individuals for now 
# for chr in {1..22}; do 
# 	gunzip ${indir}/snp_locations-chr${chr}.gz
# 	n_ind=$(cat $subj | wc -l)
# 	#for ind in $(seq 1 $n_ind); do
# 	for ind in $(seq 1 10); do
# 		echo "Working on chromosome ${chr} individual ${ind}..."
# 		prefix=${outdir}/whole-outputs/gtex-${ind}-chr${chr}
# 		start=$((${ind}*2-1))
# 		end=$((${ind}*2))
# 		tg_start=$(($n_ind*2+1))
# 		tg_end=$(($(head -1 ${indir}/classes-chr${chr} | wc -w)))
# 		cut -d' ' -f${start}-${end},${tg_start}-${tg_end} ${indir}/classes-chr${chr} > ${indir}/classes-${chr}-${ind}

# 		alleles_in=${indir}/alleles-chr${chr}
# 		cut -c${start}-${end},${tg_start}-${tg_end} ${alleles_in} > ${indir}/alleles-${ind}-chr${chr}
		
# 		classes=${indir}/classes-${chr}-${ind}
# 		snp_locations=${indir}/snp_locations-chr${chr}
# 		alleles=${indir}/alleles-${ind}-chr${chr}

# 		taskset -c 24-47 python ~/RFMix_v1.5.4/RunRFMix.py \
# 			PopPhased \
# 			${alleles} \
# 			${classes} \
# 			${snp_locations} \
# 			'~/RFMix_v1.5.4' \
# 			--forward-backward \
# 			--skip-check-input-format \
# 			-o ${prefix} &

# 		running=`ps -ef | grep "nicolerg" | grep "RunRFMix.py" | wc -l`
# 		while [ $running -gt 24 ]; do
# 			sleep 120
# 			running=`ps -ef | grep "nicolerg" | grep "RunRFMix.py" | wc -l`
# 		done 
# 	done
# done

# # clean up the outputs ================================================================================================

# # rm ${outdir}/alleles-*-chr*
# # gzip ${outdir}/snp_locations*
# # rm ${outdir}/classes-*
# # rm ${outdir}/whole-outputs/*log.txt

# cd ${outdir}/whole-outputs

# # compress whole outputs
# for chr in {1..22}; do 
# 	paste -d'\0' gtex-{1..10}-chr${chr}.0.Viterbi.txt > gtex-chr${chr}.0.Viterbi.txt &
# 	paste -d'\0' gtex-{1..10}-chr${chr}.allelesRephased0.txt > gtex-chr${chr}.allelesRephased0.txt &
# 	paste -d'\0' gtex-{1..10}-chr${chr}.0.ForwardBackward.txt > gtex-chr${chr}.0.ForwardBackward.txt &
# done
# mkdir -p merged
# mv gtex-chr* merged
# gzip merged/*

# gzip ../../rfmix_inputs/snp_locations-chr*

# generate 2 bed files for each individual (1 per haplotype) (need to make genotype track plot) ================================================

cd ${outdir}
mkdir -p bed
indir=whole-outputs/merged
ids=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/gtex-admixed0.9.txt

head ${ids} > tmp_ids.txt # I only processed the first 10 individuals in the list 

# while read ind; do

# 	if [[ -f "bed/${ind}_A.bed" ]]; then 
# 		echo '$ind has already been processed'
# 		continue 
# 	fi # skip IDs that have already been processed

# 	echo "Working on ${ind}..."
# 	echo "Collapsing results into .bed files"
# 	# generate _A.bed and _B.bed for every individual

# 	# this uses an adjusted version of one of Alicia's scripts to work on a per-chromosome basis 
# 	taskset -c 24-47 python ~/gtex-admix/rfmix_to_ancestry_bed/finalized/collapse_ancestry_single_chr.py \
# 	--rfmix ${indir}/gtex-chr19.0.Viterbi.txt.gz \
# 	--snp_locations ../rfmix_inputs/snp_locations-chr19.gz \
# 	--snp_map ../rfmix_inputs/snp-gd-bed-chr19.gz \
# 	--ind ${ind} \
# 	--ind_info ${ids} \
# 	--pop_labels 'ASN,EUR,AFR' \
# 	--out bed/${ind} \
# 	--fbk ${indir}/gtex-chr19.0.ForwardBackward.txt.gz \
# 	--fbk_threshold 0.9 &

# 	running=`ps -ef | grep "nicolerg" | grep "collapse_ancestry" | wc -l`
# 	while [ $running -gt 24 ]; do
# 		sleep 30
# 		running=`ps -ef | grep "nicolerg" | grep "collapse_ancestry" | wc -l`
# 	done

# done < tmp_ids.txt
# #done < ${ids}

# wait

# # plot karyograms
# while read ind; do
# 	python ~/gtex-admix/scripts/ancestry_pipeline/plot_karyogram.py \
# 	--bed_a bed/${ind}_A.bed \
# 	--bed_b bed/${ind}_B.bed \
# 	--ind ${ind} \
# 	--centromeres ~/gtex-admix/scripts/ancestry_pipeline//centromeres_hg19.bed \
# 	--pop_order 'ASN,EUR,AFR' \
# 	--out bed/${ind}-karyo.png &
# #done < ${ids}
# done < tmp_ids.txt

# while read ind; do
# 	python ~/gtex-admix/scripts/ancestry_pipeline/plot_karyogram.py \
# 	--bed_a bed/old/${ind}_A.bed \
# 	--bed_b bed/old/${ind}_B.bed \
# 	--ind ${ind} \
# 	--centromeres ~/gtex-admix/scripts/ancestry_pipeline//centromeres_hg19.bed \
# 	--pop_order 'ASN,EUR,AFR' \
# 	--out bed/old/${ind}-old-karyo.png &
# #done < ${ids}
# done < tmp_ids.txt

# beddir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_outputs/bed/old

# ls $beddir | grep "GTEX" | sed "s/_.*//" | sort | uniq > $beddir/tmp
# sed "s/$/_A.bed/" $beddir/tmp > $beddir/tmp1
# sed "s/$/_B.bed/" $beddir/tmp > $beddir/tmp2
# paste $beddir/tmp1 $beddir/tmp2 > $beddir/gtex.bed.list.txt
# rm $beddir/tmp*

# # calculate global ancestry
# echo "Calculating global ancestry"
# python /users/nicolerg/gtex-admix/scripts/ancestry_pipeline/lai_global.py \
# --bed_list ${beddir}/gtex.bed.list.txt \
# --ind_list tmp_ids.txt \
# --pops 'ASN,EUR,AFR,UNK' \
# --out ${beddir}/gtex_ai_global.txt \
# --indir ${beddir}

# mkdir -p $beddir/master

beddir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_outputs/bed/old

# # combine bed files into master ancestry files ================================================================================================
# for chr in {1..22}; do
# 	outfile=${beddir}/ancestry.chr${chr}.bed
# 	for file in `ls ${beddir} | grep "GTEX-"`; do
# 		ind=`echo ${file} | sed "s/\..*//"`
# 		grep -E "^${chr}	" ${beddir}/${file} | cut -f 1,2,3,4 | sed "s/$/	${ind}/" | sed "s/^/chr/" >> $outfile
# 	done
# done

# make master ancestry files from haplotype bed files
# this takes a while to run
for chr in {1..22}; do
	taskset -c 24-47 python ${srcdir}/interpolate-local-anc.py ${chr} &

	running=`ps -ef | grep "nicolerg" | grep "interpolate-local-anc" | wc -l`
	while [ $running -gt 5 ]; do
		sleep 30
		running=`ps -ef | grep "nicolerg" | grep "interpolate-local-anc" | wc -l`
	done
done


