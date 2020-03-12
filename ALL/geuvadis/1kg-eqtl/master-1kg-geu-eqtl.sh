# GEUVADIS eQTL discovery with GTEx pipeline

indir="/mnt/lab_data/montgomery/nicolerg/geuvadis"
outdir=${indir}

exprdir=${indir}/expression
peerdir=${indir}/peer 
covdir=${indir}/cov
pcadir=${indir}/pca

vcfin="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/1000genomes-phased-hg19"

# first, whole-genome things

annotation_gtf="${exprdir}/gencode.v19.genes.v7.patched_contigs.gtf" # custom GTEx GTF file compatible with normalize_expression.py

# # genotype PCA ---------------------------------------------------------------------------

# perform genotype PCA for EUR and AFR separately
if [ ! -e "${pcadir}/tg.EUR.subset.recode.vcf.gz" ]; then 
	grep "AFR" ${covdir}/geuvadis-afr-cov.tsv | cut -f1 > ${vcfdir}/afr-subjs.txt
	grep "EUR" ${covdir}/geuvadis-eur-cov.tsv | cut -f1 > ${vcfdir}/eur-subjs.txt
	sed -i '/NA19118/d' ${vcfdir}/afr-subjs.txt # remove outlier
	vcftools --gzvcf "/mnt/lab_data/montgomery/nicolerg/pca-subset/tg.master.subset.vcf.gz" --out ${pcadir}/tg.EUR.subset --keep ${vcfdir}/eur-subjs.txt --recode --recode-INFO-all &
	vcftools --gzvcf "/mnt/lab_data/montgomery/nicolerg/pca-subset/tg.master.subset.vcf.gz" --out ${pcadir}/tg.AFR.subset --keep ${vcfdir}/afr-subjs.txt --recode --recode-INFO-all &
	wait
	bgzip -f ${pcadir}/tg.EUR.subset.recode.vcf
	bgzip -f ${pcadir}/tg.AFR.subset.recode.vcf
fi

# normalize expression ---------------------------------------------------------------------------

expr=${exprdir}/GD660.GeneQuantRPKM.filt.gct.gz

# The expression data are normalized as follows: 
# (i) expression values are quantile normalized to the average empirical distribution observed across samples; 
# (ii) for each gene, expression values are inverse quantile normalized to a standard normal distribution across samples.
# Using these settings, genes are selected based on expression thresholds of >0 RPKM in ≥90% of samples

if [ ! -e "${exprdir}/GV-EUR.expression.bed.gz" -a ! -e "${exprdir}/GV-EUR.expression.bed" ]; then
	#conda create -n py35 python=3.5
	source activate py35

	prefix="GV-all"
	nsamp=`zcat ${expr} | head -2 | tail -1 | sed "s/.*	//"`
	minsamp=$((nsamp * 9 / 10))

	taskset -c 30-47 python /users/nicolerg/gtex-admix/scripts/LAVA/1kg-eqtl/normalize_expression_rpkm.py \
	${expr} \
	${annotation_gtf} \
	'/mnt/lab_data/montgomery/shared/1KG/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' \
	${prefix} \
	-o ${exprdir} \
	--expression_threshold 0.0 \
	--min_samples ${minsamp} 

	source deactivate

	# now split them into EUR and AFR

	met='/users/nicolerg/gtex-admix/metadata/geuvadis-cov.tsv'
	for pop in "EUR" "AFR"; do
		zcat ${exprdir}/GV-all.expression.bed.gz | cut -f1,2,3,4 > ${exprdir}/GV-${pop}.expression.bed
		for sub in `grep "$pop" $met | cut -f1`; do 
			# exclude outlier
			if [ "$sub" == "NA19118" ]; then
				continue
			fi
			col=`zcat ${exprdir}/GV-all.expression.bed.gz | head -1 | sed 's/	/\n/g' | grep -n "$sub" | sed "s/:.*//"`
			if [ "$col" != "" ]; then
				zcat ${exprdir}/GV-all.expression.bed.gz | cut -f$col > ${exprdir}/tmp
				paste ${exprdir}/GV-${pop}.expression.bed ${exprdir}/tmp > ${exprdir}/tmp.2
				rm ${exprdir}/GV-${pop}.expression.bed 
				mv ${exprdir}/tmp.2 ${exprdir}/GV-${pop}.expression.bed
			fi
		done
	done

fi

# # to instead normalize the populations separately:
# source activate py35

# for pop in "EUR" "AFR"
# do
# 	norm="${exprdir}/GEUV-${pop}-expr.gct.gz"
# 	prefix="GV-${pop}"
# 	nsamp=`zcat ${norm} | head -2 | tail -1 | sed "s/.*	//"`
# 	minsamp=$((nsamp / 2))

# 	taskset -c 30-47 python /users/nicolerg/gtex-admix/scripts/LAVA/1kg-eqtl/normalize_expression_rpkm.py \
# 	${norm} \
# 	${annotation_gtf} \
# 	${vcfdir}/tg.${pop}.chr${chr}.recode.vcf.gz \
# 	${prefix} \
# 	-o ${exprdir} \
# 	--expression_threshold 1.0 \
# 	--min_samples ${minsamp} &
# 	# This step will generate 5 files:
# 	# ${prefix}.expression.bed.gz
# 	# ${prefix}.expression.bed.gz.tbi
# 	# ${prefix}.expression.fpkm.bed.gz
# 	# ${prefix}.expression.fpkm.bed.gz.tbi
# 	# ${prefix}.expression.txt
# done

# source deactivate

# wait

# # IF K (N PEER factors) HAS NOT YET BEEN SELECTED ---------------------------------------------------------------------------

# # select optimal number of PEER factors
# for pop in "EUR" "AFR"
# do
# 	#run chr22 with 60 PEER factors
# 	prefix="GV-${pop}"

# 	num_peer=60
# 	taskset -c 30-47 Rscript /users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/run_PEER.R ${exprdir}/${prefix}.expression.txt ${prefix} ${num_peer} -o ${peerdir} &
# 	# This will generate 3 files:
# 	# ${prefix}_PEER_residuals.txt
# 	# ${prefix}_PEER_alpha.txt
# 	# ${prefix}_PEER_covariates.txt
# done

# wait

# # make covariates file
# Rscript /users/nicolerg/gtex-admix/scripts/geuvadis/1kg-eqtl/genotype-pca.R

# # add in PEER factors to see how it affects discovery rate
# for pop in "EUR" "AFR"
# do
# 	basecov="${covdir}/${pop}.covariates.tsv"
# 	if [ "$pop" == "EUR" ]
# 	then
# 		base=4
# 		vcf=$eur_vcf
# 		baseexpr=${exprdir}/GV-EUR.expression.bed.gz
# 	else
# 		base=5
# 		vcf=$afr_vcf
# 		baseexpr=${exprdir}/GV-AFR.expression.bed.gz
# 	fi

# 	# select only expression in this chromosome
# 	tbi="${exprdir}/GV-${pop}.chr${chr}.expr.bed.gz.tbi"
# 	if [ ! -e "$tbi" ]; then
# 		zcat ${baseexpr} | head -1 > "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
# 		zgrep -E "^chr${chr}	" ${baseexpr} >> "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
# 		bgzip "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
# 		tabix -p bed "${exprdir}/GV-${pop}.chr${chr}.expr.bed.gz"
# 	fi
# 	expr="${exprdir}/GV-${pop}.chr${chr}.expr.bed.gz"

# 	for K in 0 1 3 5 7 10 13 15 20
# 	do
# 		i=$(($K+$base))
# 		sed -n "1,${i}p" ${basecov} > ${covdir}/${pop}.cov.${K}.peer.txt

# 		for j in $(seq 1 30); do
# 			taskset -c 30-47 /users/nicolerg/FastQTL/bin/fastQTL \
# 			--vcf ${vcf} \
# 			--bed ${expr} \
# 			--out ${eqtldir}/${pop}.chr${chr}.nominals.${K}.chunk${j}.txt.gz \
# 			--window 1e6 \
# 			--cov ${covdir}/${pop}.cov.${K}.peer.txt \
# 			--chunk $j 30 &
# 		done

# 	done

# 	wait

# 	for K in 0 1 3 5 7 10 13 15 20; do 

# 		rm ${covdir}/${pop}.cov.${K}.peer.txt

# 		# merge chunks
# 		for file in ${eqtldir}/${pop}.chr${chr}.nominals.${K}.chunk{1..30}.txt.gz
# 		do 
# 			zcat $file >> ${eqtldir}/${pop}.chr${chr}.nominals.${K}.merged.txt
# 		done
# 		rm ${eqtldir}/${pop}.chr${chr}.nominals.${K}.chunk*
# 		gzip ${eqtldir}/${pop}.chr${chr}.nominals.${K}.merged.txt

# 	done

# done

# # get # eQTLs in each iteration
# for pop in "EUR" "AFR"; do
# 	echo "n_peer"	"n_eQTL_genes" > ${peerdir}/${pop}.chr${chr}.n.eqtl.txt
# 	for K in 0 1 3 5 7 10 13 15 20; do
# 		file=${eqtldir}/${pop}.chr${chr}.nominals.${K}.merged.txt
# 		# correct for multiple testing
# 		ntests=`zcat $file | wc -l`
# 		p=`echo $(echo "0.05/$ntests" | bc -l )`
# 		#zcat ${eqtldir}/${pop}.chr${chr}.nominals.${K}.merged.txt | awk -F' ' -v p="$p" '{if($4<=p){print $1}}' | sort | uniq | wc -l
# 		nsig=`zcat $file | awk -F' ' -v p="$p" '{if($4<=p){print $1}}' | sort | uniq | wc -l`
# 		echo $K	$nsig >> ${peerdir}/${pop}.chr${chr}.n.eqtl.txt
# 	done
# done

# look at counts to pick K for PEER ---> use K=13

# # ------------------------------------------------------------------------------------------------------------------------------------------------------

# make covariates file with selected number of PEER covariates

# 15 PEER factors for AFR
# 60 PEER factors for EUR

if [ ! -e "${covdir}/EUR.full.covariates.tsv" -o ! -e "${covdir}/AFR.full.covariates.tsv" ]; then # if final covariates file doesn't exist yet 
	# first for EUR 
	K=60
	prefix="GV-EUR"
	taskset -c 30-47 Rscript /users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/run_PEER.R ${exprdir}/${prefix}.expression.bed ${prefix} ${K} -o ${peerdir} &
	# for AFR
	K=15
	prefix="GV-AFR"
	taskset -c 30-47 Rscript /users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/run_PEER.R ${exprdir}/${prefix}.expression.bed ${prefix} ${K} -o ${peerdir} &

	wait

	# make covariates file
	Rscript /users/nicolerg/gtex-admix/scripts/geuvadis/1kg-eqtl/genotype-pca.R

	# # this only applies if you calculate more PEER factors than necessary
	# # for EUR 
	# basecov=${covdir}/EUR.covariates.tsv
	# l=`grep -n "InferredCov60" ${basecov} | sed "s/\:.*//"`
	# head -$l ${basecov} > ${covdir}/EUR.full.covariates.tsv
	# # for AFR
	# basecov=${covdir}/AFR.covariates.tsv
	# l=`grep -n "InferredCov15" ${basecov} | sed "s/\:.*//"`
	# head -$l ${basecov} > ${covdir}/AFR.full.covariates.tsv

	cp ${covdir}/EUR.covariates.tsv ${covdir}/EUR.full.covariates.tsv
	cp ${covdir}/AFR.covariates.tsv ${covdir}/AFR.full.covariates.tsv

	gzip ${exprdir}/GV-EUR.expression.bed
	gzip ${exprdir}/GV-AFR.expression.bed
fi

# then per-chromosome things

# # EUR into AFR =====================================================================================================================

# eqtldir=${indir}/eqtl/EUR-into-AFR
# vcfdir=${indir}/genotype/EUR-into-AFR

# # make filtered VCFs
# if [ ! -e "${vcfdir}/tg.EUR.chr22.recode.vcf.gz" ]; then
# 	for chr in {1..22}; do 
# 		bash '/users/nicolerg/gtex-admix/scripts/LAVA/1kg-eqtl/make-vcfs.sh' ${chr} 'EUR' &
# 	done
# fi

# wait

# for chr in {1..22}; do

# 	eur_vcf=${vcfdir}/tg.EUR.chr${chr}.recode.vcf.gz
# 	afr_vcf=${vcfdir}/tg.AFR.chr${chr}.recode.vcf.gz

# 	# run FASTQTL

# 	nchunks=50
# 	for pop in "AFR" "EUR"; do
# 		cov="${covdir}/${pop}.full.covariates.tsv"

# 		# select only expression in this chromosome
# 		baseexpr=${exprdir}/GV-${pop}.expression.bed.gz
# 		expr="${exprdir}/GV-${pop}.chr${chr}.expr.bed.gz"
# 		if [ ! -e "$expr" ]; then
# 			zcat "${baseexpr}" | head -1 > "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
# 			#echo "here"
# 			zgrep -E "^chr${chr}	" ${baseexpr} >> "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
# 			bgzip "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
# 			tabix -p bed "$expr"
# 		fi

# 		if [ "$pop" == "EUR" ]; then
# 			vcf=$eur_vcf	
# 		else
# 			vcf=$afr_vcf
# 		fi

# 		# nominal pass
# 		for j in $(seq 1 $nchunks); do
# 			taskset -c 30-47 /users/nicolerg/FastQTL/bin/fastQTL \
# 			--vcf ${vcf} \
# 			--bed ${expr} \
# 			--window 1e6 \
# 			--out ${eqtldir}/${pop}.chr${chr}.chunk${j}.txt.gz \
# 			--cov ${cov} \
# 			--chunk $j $nchunks &

# 			# limit the number of jobs submitted at once 
# 			running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
# 			while [ $running -gt 30 ] # run max 30 jobs at a time
# 			do
# 				sleep 10
# 				running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
# 			done

# 		done

# 		# permutation pass
# 		for j in $(seq 1 $nchunks); do
# 			taskset -c 30-47 /users/nicolerg/FastQTL/bin/fastQTL \
# 			--vcf ${vcf} \
# 			--bed ${expr} \
# 			--window 1e6 \
# 			--out ${eqtldir}/${pop}.chr${chr}.permuted.chunk${j}.txt.gz \
# 			--cov ${cov} \
# 			--permute 1000 10000 \
# 			--chunk $j $nchunks &

# 			# limit the number of jobs submitted at once 

# 			running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
# 			while [ $running -gt 12 ] 
# 			do
# 				sleep 10
# 				running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
# 			done

# 		done

# 	done

# 	wait

# 	# merge chunks - nominal pass
# 	for pop in "AFR" "EUR"; do
# 		if [ ! -e "${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt.gz" ]; then
# 			echo "gene_id variant tss_distance pval_nominal slope" > ${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt
# 			for j in $(seq 1 $nchunks); do
# 				file=${eqtldir}/${pop}.chr${chr}.chunk${j}.txt.gz
# 				zcat $file >> ${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt
# 				rm $file 
# 			done
# 			gzip ${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt
# 		fi
# 	done

# 	# merge chunks - permutation pass
# 	for pop in "AFR" "EUR"; do
# 		if [ ! -e "${eqtldir}/${pop}.chr${chr}.geu.permuted.txt.gz" ]; then
# 			echo "gene_id num_var shape1_mle shape2_mle dummy variant tss_distance pval_nominal slope pval_perm pval_beta" > ${eqtldir}/${pop}.chr${chr}.geu.permuted.txt
# 			for j in $(seq 1 $nchunks); do
# 				file=${eqtldir}/${pop}.chr${chr}.permuted.chunk${j}.txt.gz
# 				zcat $file >> ${eqtldir}/${pop}.chr${chr}.geu.permuted.txt
# 				rm $file 
# 			done
# 			gzip ${eqtldir}/${pop}.chr${chr}.geu.permuted.txt
# 		fi
# 	done

# done

# # merge chr into whole genome - nominal pass
# for pop in "EUR" "AFR"; do
# 	if [ ! -e "${eqtldir}/${pop}.geu.allpairs.txt.gz" ]; then
# 		zcat ${eqtldir}/${pop}.chr1.geu.allpairs.txt.gz | head -1 > ${eqtldir}/${pop}.geu.allpairs.txt
# 		for file in ${eqtldir}/${pop}.chr{1..22}.geu.allpairs.txt.gz; do 
# 			zcat $file | sed -e '1d' >> ${eqtldir}/${pop}.geu.allpairs.txt
# 			# rm $file
# 		done 
# 		gzip ${eqtldir}/${pop}.geu.allpairs.txt
# 	fi
# done

# # merge chr into whole genome - permutation pass
# for pop in "EUR" "AFR"; do
# 	if [ ! -e "${eqtldir}/${pop}.geu.permuted.txt.gz" ]; then
# 		zcat ${eqtldir}/${pop}.chr1.geu.permuted.txt.gz | head -1 > ${eqtldir}/${pop}.geu.permuted.txt
# 		for file in ${eqtldir}/${pop}.chr{1..22}.geu.permuted.txt.gz; do 
# 			zcat $file | sed -e '1d' >> ${eqtldir}/${pop}.geu.permuted.txt
# 			# rm $file
# 		done 
# 		gzip ${eqtldir}/${pop}.geu.permuted.txt
# 	fi
# done

# # AFR into EUR =====================================================================================================================

eqtldir=${indir}/eqtl/AFR-into-EUR
vcfdir=${indir}/genotype/AFR-into-EUR

# make filtered VCFs
if [ ! -e "${vcfdir}/tg.EUR.chr22.recode.vcf.gz" ]; then
	for chr in {1..22}; do 
		bash '/users/nicolerg/gtex-admix/scripts/LAVA/1kg-eqtl/make-vcfs.sh' ${chr} 'AFR' &
	done
fi

wait

for chr in {1..22}; do

	eur_vcf=${vcfdir}/tg.EUR.chr${chr}.recode.vcf.gz
	afr_vcf=${vcfdir}/tg.AFR.chr${chr}.recode.vcf.gz

	# run FASTQTL

	nchunks=50
	for pop in "AFR" "EUR"; do
		cov="${covdir}/${pop}.full.covariates.tsv"

		# select only expression in this chromosome
		baseexpr=${exprdir}/GV-${pop}.expression.bed.gz
		expr="${exprdir}/GV-${pop}.chr${chr}.expr.bed.gz"
		if [ ! -e "$expr" ]; then
			zcat "${baseexpr}" | head -1 > "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
			#echo "here"
			zgrep -E "^chr${chr}	" ${baseexpr} >> "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
			bgzip "${exprdir}/GV-${pop}.chr${chr}.expr.bed"
			tabix -p bed "$expr"
		fi

		if [ "$pop" == "EUR" ]; then
			vcf=$eur_vcf	
		else
			vcf=$afr_vcf
		fi

		# nominal pass
		for j in $(seq 1 $nchunks); do
			taskset -c 30-47 /users/nicolerg/FastQTL/bin/fastQTL \
			--vcf ${vcf} \
			--bed ${expr} \
			--window 1e6 \
			--out ${eqtldir}/${pop}.chr${chr}.chunk${j}.txt.gz \
			--cov ${cov} \
			--chunk $j $nchunks &

			# limit the number of jobs submitted at once 
			running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
			while [ $running -gt 30 ] # run max 30 jobs at a time
			do
				sleep 10
				running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
			done

		done

		# permutation pass
		for j in $(seq 1 $nchunks); do
			taskset -c 30-47 /users/nicolerg/FastQTL/bin/fastQTL \
			--vcf ${vcf} \
			--bed ${expr} \
			--window 1e6 \
			--out ${eqtldir}/${pop}.chr${chr}.permuted.chunk${j}.txt.gz \
			--cov ${cov} \
			--permute 1000 10000 \
			--chunk $j $nchunks &

			# limit the number of jobs submitted at once 

			running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
			while [ $running -gt 12 ] 
			do
				sleep 10
				running=`ps -ef | grep "nicolerg" | grep "fastQTL" | wc -l`
			done

		done

	done

	wait

	# merge chunks - nominal pass
	for pop in "AFR" "EUR"; do
		if [ ! -e "${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt.gz" ]; then
			echo "gene_id variant tss_distance pval_nominal slope" > ${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt
			for j in $(seq 1 $nchunks); do
				file=${eqtldir}/${pop}.chr${chr}.chunk${j}.txt.gz
				zcat $file >> ${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt
				rm $file 
			done
			gzip ${eqtldir}/${pop}.chr${chr}.geu.allpairs.txt
		fi
	done

	# merge chunks - permutation pass
	for pop in "AFR" "EUR"; do
		if [ ! -e "${eqtldir}/${pop}.chr${chr}.geu.permuted.txt.gz" ]; then
			echo "gene_id num_var shape1_mle shape2_mle dummy variant tss_distance pval_nominal slope pval_perm pval_beta" > ${eqtldir}/${pop}.chr${chr}.geu.permuted.txt
			for j in $(seq 1 $nchunks); do
				file=${eqtldir}/${pop}.chr${chr}.permuted.chunk${j}.txt.gz
				zcat $file >> ${eqtldir}/${pop}.chr${chr}.geu.permuted.txt
				rm $file 
			done
			gzip ${eqtldir}/${pop}.chr${chr}.geu.permuted.txt
		fi
	done

done

# merge chr into whole genome - nominal pass
for pop in "EUR" "AFR"; do
	if [ ! -e "${eqtldir}/${pop}.geu.allpairs.txt.gz" ]; then
		zcat ${eqtldir}/${pop}.chr1.geu.allpairs.txt.gz | head -1 > ${eqtldir}/${pop}.geu.allpairs.txt
		for file in ${eqtldir}/${pop}.chr{1..22}.geu.allpairs.txt.gz; do 
			zcat $file | sed -e '1d' >> ${eqtldir}/${pop}.geu.allpairs.txt
			# rm $file
		done 
		gzip ${eqtldir}/${pop}.geu.allpairs.txt
	fi
done

# merge chr into whole genome - permutation pass
for pop in "EUR" "AFR"; do
	if [ ! -e "${eqtldir}/${pop}.geu.permuted.txt.gz" ]; then
		zcat ${eqtldir}/${pop}.chr1.geu.permuted.txt.gz | head -1 > ${eqtldir}/${pop}.geu.permuted.txt
		for file in ${eqtldir}/${pop}.chr{1..22}.geu.permuted.txt.gz; do 
			zcat $file | sed -e '1d' >> ${eqtldir}/${pop}.geu.permuted.txt
			# rm $file
		done 
		gzip ${eqtldir}/${pop}.geu.permuted.txt
	fi
done

# now we have two sets of results: 
# 1. /mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR-into-AFR/[AFR,EUR].geu.[allpairs,permuted].txt.gz
# 2. /mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR/[AFR,EUR].geu.[allpairs,permuted].txt.gz

ID eQTLs (qval of pval_beta + 5% FDR cutoff); plot EUR vs AFR effect sizes
Rscript /users/nicolerg/gtex-admix/scripts/geuvadis/geu-custom-effect-sizes.R '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR-into-AFR/'
Rscript /users/nicolerg/gtex-admix/scripts/geuvadis/geu-custom-effect-sizes.R '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR/'

