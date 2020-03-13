#!/bin/bash
# Nicole Gay

# updated 25 April 2019
# local ancestry pipeline
# 1. format inputs for local ancestry imputation (RFMix)
# 2. run RFMix
# 3. generate BED files
# 4. calculate global ancestry 
# 5. make master ancestry files from haplotype bed files 
# 6. make local covariate files 

set -e

srcdir=~/gtex-admix/scripts/rfmix_to_ancestry_bed # path to scripts 
indir=/mnt/lab_data/montgomery/nicolerg/genotypes-hg19 # path to VCF files (all hg19 for compatibility with genetic map)
outdir=/mnt/lab_data/montgomery/nicolerg/rfmix-all # output directory
gpheno=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt # admixed subject phenotype file 
tdict=/users/nicolerg/gtex-admix/metadata/1kg-subject-popcodes.tsv # TSV file that matches 1000 Genomes subject to subpopulations
subj=/mnt/lab_data/montgomery/nicolerg/rfmix-all/gtex-admix0.9.txt # list of GTEx/admixed IDs, one ID per line (must match column names in admixed VCF)

# format inputs ==========================================================================================

for chr in {1..22}; do 
	tg=${indir}/1000genomes-phased-hg19/tg.chr${chr}.vcf.gz # 1000 genomes VCF, hg19
	gtex=${indir}/gtex-phased-MAF0-hg19/gtex.chr${chr}.vcf.gz # GTEx VCF, hg19

	echo "Generating VCF for chr${chr}..."
	Rscript ${srcdir}/format-rfmix.R ${subj} ${tdict} ${gtex} ${tg} ${outdir} ${chr}
	gzip ${outdir}/alleles-chr${chr}
	gzip ${outdir}/snp-bed-chr${chr}
done

# condense class and subject ID files
mv ${outdir}/$(ls ${outdir} | grep "classes-chr" | head -1) ${outdir}/classes
rm -f ${outdir}/classes-chr*
mv ${outdir}/$(ls ${outdir} | grep "subject-ids-chr" | head -1) ${outdir}/subject-ids
# mv ${outdir}/subject-ids-chr1 ${outdir}/subject-ids
rm -f ${outdir}/subject-ids-chr*

# convert bp to genetic distance ============================================================================

echo "Converting bp to cM..."
for chr in {1..22}; do
	snp_bed="${outdir}/snp-bed-chr${chr}.gz" # output from previous step
	mapfile="/mnt/lab_data/montgomery/nicolerg/RESOURCES/genetic-map/hapmap_genetic_maps/genetic_map_GRCh37_chr${chr}.txt.gz" # genetic map, required to convert bp to cM (hg19)

	python ${srcdir}/bp-to-cM.py ${snp_bed} ${mapfile} ${outdir}/snp-gd-bed-chr${chr}.gz ${outdir}/snp_locations-chr${chr}.gz &

	# limit the number of jobs running in parallel
	running=`ps -ef | grep "nicolerg" | grep "bp-to-cM" | wc -l`
	while [ $running -gt 8 ]; do
		sleep 30
		running=`ps -ef | grep "nicolerg" | grep "bp-to-cM" | wc -l`
	done 
done

rm -f ${outdir}/snp-bed-chr*.gz

wait

# run RFMix ================================================================================================

mkdir -p $outdir

# this block runs the algorithm separately for each admixed individual
# it iteratively parses the input files to include only one admixed individual and all reference panels
for chr in {1..22}; do 
	gunzip ${outdir}/snp_locations-chr${chr}.gz
	n_ind=$(cat $subj | wc -l)
	for ind in $(seq 1 $n_ind); do
		echo "Working on chromosome ${chr} individual ${ind}..."
		prefix=${outdir}/whole-outputs/gtex-${ind}-chr${chr}
		start=$((${ind}*2-1))
		end=$((${ind}*2))
		tg_start=$(($n_ind*2+1))
		tg_end=$(($(head -1 ${outdir}/classes | wc -w)))
		cut -d' ' -f${start}-${end},${tg_start}-${tg_end} ${outdir}/classes > ${outdir}/classes-${ind}

		alleles_in=${outdir}/alleles-chr${chr}.gz
		# zcat ${alleles_in} | cut -c${start}-${end},53-1474 > ${outdir}/alleles-${ind}-chr${chr}
		zcat ${alleles_in} | cut -c${start}-${end},${tg_start}-${tg_end} > ${outdir}/alleles-${ind}-chr${chr}
		
		classes=${outdir}/classes-${ind}
		snp_locations=${outdir}/snp_locations-chr${chr}
		alleles=${outdir}/alleles-${ind}-chr${chr}

		taskset -c 24-47 python ~/RFMix_v1.5.4/RunRFMix.py \
			PopPhased \
			${alleles} \
			${classes} \
			${snp_locations} \
			'~/RFMix_v1.5.4' \
			--forward-backward \
			--skip-check-input-format \
			-o ${prefix} &
		
		running=`ps -ef | grep "nicolerg" | grep "RunRFMix.py" | wc -l`
		while [ $running -gt 25 ]; do
			sleep 120
			running=`ps -ef | grep "nicolerg" | grep "RunRFMix.py" | wc -l`
		done 
	done
done

wait

# clean up the outputs ================================================================================================

rm ${outdir}/alleles-*-chr*
gzip ${outdir}/snp_locations*
rm ${outdir}/classes-*
# rm ${outdir}/whole-outputs/*log.txt

# compress whole outputs
for chr in {1..22}; do 
	paste -d'\0' ${outdir}/whole-outputs/gtex-*-chr${chr}.0.Viterbi.txt > ${outdir}/gtex-chr${chr}.0.Viterbi.txt ### THIS COLLAPSES SUBJECTS IN THE WRONG ORDER
	paste -d'\0' ${outdir}/whole-outputs/gtex-*-chr${chr}.allelesRephased0.txt > ${outdir}/gtex-chr${chr}.allelesRephased0.txt
	paste -d'\0' ${outdir}/whole-outputs/gtex-*-chr${chr}.0.ForwardBackward.txt > ${outdir}/gtex-chr${chr}.0.ForwardBackward.txt

	gzip ${outdir}/whole-outputs/gtex-chr${chr}.0.Viterbi.txt
	gzip ${outdir}/whole-outputs/gtex-chr${chr}.allelesRephased0.txt
	gzip ${outdir}/whole-outputs/gtex-chr${chr}.0.ForwardBackward.txt
done

# mkdir -p ${outdir}/inputs
# mv ${outdir}/classes ${outdir}/subject-ids ${outdir}/alleles-* ${outdir}/snp_locations-* ${outdir}/inputs 

# generate 2 bed files for each individual (1 per haplotype) (need to make genotype track plot) ================================================

beddir=${outdir}/bed
mkdir -p $beddir
indir=${outdir}/whole-outputs
ids=${outdir}/remaining_ids.txt

while read ind
do
	if [[ -f "${beddir}/${ind}_A.bed" ]]; then 
		echo '$ind has already been processed'
		continue 
	fi # skip IDs that have already been processed

	echo "Working on ${ind}..."
	echo "Collapsing results into .bed files"
	# generate _A.bed and _B.bed for every individual

	# this uses an adjusted version of one of Alicia's scripts to work on a per-chromosome basis 
	taskset -c 24-47 python ${srcdir}/collapse_ancestry_single_chr.py \
	--rfmix ${indir}/gtex-chr19.0.Viterbi.txt.gz \
	--snp_locations ${outdir}/inputs/snp_locations-chr19.gz \
	--snp_map /mnt/lab_data/montgomery/nicolerg/rfmix-all/snp-gd-bed-chr19.gz \
	--chr 'chr19' \
	--ind ${ind} \
	--ind_info ${outdir}/inputs/subject-ids \
	--pop_labels 'ASN,EUR,AFR' \
	--out ${beddir}/${ind} \
	--fbk ${indir}/gtex-chr19.0.ForwardBackward.txt.gz \
	--fbk_threshold 0.9 &

	running=`ps -ef | grep "nicolerg" | grep "collapse_ancestry" | wc -l`
	while [ $running -gt 24 ]; do
		sleep 30
		running=`ps -ef | grep "nicolerg" | grep "collapse_ancestry" | wc -l`
	done

done < ${ids}

wait

# while read ind
# do
# 	echo "Plotting karyotypes"
# 	# plot karyotypes
# 	python ${srcdir}/plot_karyogram.py \
# 	--bed_a "${outdir}/${ind}_A.bed" \
# 	--bed_b "${outdir}/${ind}_B.bed" \
# 	--ind ${ind} \
# 	--centromeres "${srcdir}/centromeres_hg19.bed" \
# 	--pop_order 'ASN,EUR,AFR' \
# 	--out "${outdir}/${ind}-karyo.png" &
# 	echo ""
# done < ${ids}

# make list of bed files: two file names per line, corresponding to haps A and B for the same ind --> gtex.bed.list.txt
# NOTE: if combining with other subjects, do it at this point: put all haplotype-level bed files in a single folder
beddir=/mnt/lab_data/montgomery/nicolerg/rfmix-all/combined-bed

ls $beddir | sed "s/_.*//" | sort | uniq > $beddir/tmp
sed "s/$/_A.bed/" $beddir/tmp > $beddir/tmp1
sed "s/$/_B.bed/" $beddir/tmp > $beddir/tmp2
paste $beddir/tmp1 $beddir/tmp2 > $beddir/gtex.bed.list.txt
rm $beddir/tmp*

# calculate global ancestry
echo "Calculating global ancestry"
python /users/nicolerg/gtex-admix/scripts/ancestry_pipeline/lai_global.py \
--bed_list ${beddir}/gtex.bed.list.txt \
--ind_list /mnt/lab_data/montgomery/nicolerg/rfmix-all/all_ids_sorted.txt \
--pops 'ASN,EUR,AFR,UNK' \
--out ${beddir}/gtex_ai_global.txt

# # plot haplotype block 
# srcdir="/users/nicolerg/gtex-admix/scripts/rfmix_analysis"
# Rscript ${srcdir}/plot-haplotypes.R

mkdir -p $beddir/master

# combine bed files into master ancestry files ================================================================================================
for chr in {1..22}; do
	outfile=${beddir}/master/ancestry.chr${chr}.bed
	for file in `ls ${beddir} | grep "GTEX-"`; do
		ind=`echo ${file} | sed "s/\..*//"`
		grep -E "^${chr}	" ${beddir}/${file} | cut -f 1,2,3,4 | sed "s/$/	${ind}/" | sed "s/^/chr/" >> $outfile
	done
done

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


