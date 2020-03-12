# GTEx v7 eQTL discovery pipeline with scripts from https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl with updated parameters from Francois

base="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8"
gctbase="/mnt/lab_data/montgomery/nicolerg/gtex-sampleannot"
annotation_gtf="${base}/references/gencode.v26.GRCh38.genes.gtf"
vcf="${gctbase}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
outdir="/mnt/lab_data/montgomery/nicolerg/eqtl"

# The expression data are normalized as follows: 
# (i) expression values are quantile normalized to the average empirical distribution observed across samples; 
# (ii) for each gene, expression values are inverse quantile normalized to a standard normal distribution across samples.
# Using these settings, genes are selected based on expression thresholds of >0.1 RPKM in ≥20 samples and >6 reads in ≥20 samples

#conda create -n py35 python=3.5
source activate py35

for sys in `ls ${gctbase} | grep "_gene_reads" | sed 's/_gene.*//'`
do
	counts="${gctbase}/${sys}_gene_reads.gct.gz"
	norm="${gctbase}/${sys}_gene_tpm.gct.gz"
	prefix="${sys}"

	nsamp=`zcat ${counts} | head -3 | tail -1 | wc -w`
	if [ "${nsamp}" -ge "70" ] # only do this if there are at least 70 samples for the tissue 
	then 
		taskset -c 12-47 python /users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/normalize_expression.py \
		${norm} \
		${counts} \
		${annotation_gtf} \
		${vcf} \
		${prefix} \
		-o ${outdir} \
		--expression_threshold 0.1 \
		--count_threshold 6 \
		--min_samples 20 &
		# This step will generate 5 files:
		# ${prefix}.expression.bed.gz
		# ${prefix}.expression.bed.gz.tbi
		# ${prefix}.expression.fpkm.bed.gz
		# ${prefix}.expression.fpkm.bed.gz.tbi
		# ${prefix}.expression.txt
	fi
done

wait

source deactivate

# ---------------------------------------------------------------------------

for prefix in `ls ${outdir} | sed 's/\..*//' | sort | uniq | sed '/other_covariates/d'`
do
	file=`ls ${gctbase} | grep ${prefix} | grep "_gene_reads"`
	nsamp=`zcat ${gctbase}/${file} | head -3 | tail -1 | wc -w`

	if [ "${nsamp}" -lt "150" ]
	then
		num_peer=15
	elif [ "${nsamp}" -ge "150" ] && [ "${nsamp}" -lt "250" ]
	then
		num_peer=30
	elif [ "${nsamp}" -ge "250" ] && [ "${nsamp}" -lt "350" ]
	then
		num_peer=45
	else
		num_peer=60
	fi

	taskset -c 12-47 Rscript /users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/run_PEER.R ${outdir}/${prefix}.expression.txt ${prefix} ${num_peer} -o ${outdir} &
	# This will generate 3 files:
	# ${prefix}_PEER_residuals.txt
	# ${prefix}_PEER_alpha.txt
	# ${prefix}_PEER_covariates.txt
done

wait

# ---------------------------------------------------------------------------

# make all_covariates files

# add sex, platform, and pcr to covariate files
# reformat covariate files
outdir="/mnt/lab_data/montgomery/nicolerg/eqtl"
output="/mnt/lab_data/montgomery/nicolerg/eqtl/other_covariates.txt"
temp1=${outdir}/temp1.txt
temp2=${outdir}/temp2.txt
temp3=${outdir}/temp3.txt

cov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates"

for file in `ls ${cov}`
do 
	if [ ! -e ${temp1} ]
	then
		head -1 ${cov}/${file} > ${temp1}
		tail -3 ${cov}/${file} >> ${temp1}
	else
		head -1 ${cov}/${file} > ${temp2}
		tail -3 ${cov}/${file} >> ${temp2}
		paste ${temp1} ${temp2} > ${temp3}
		rm ${temp1}
		mv ${temp3} ${temp1}
	fi
done

mv ${temp1} ${output}
rm ${temp2}

for file in `ls ${outdir} | grep "PEER_covariates"`
do
    prefix=`echo $file | sed 's/\..*//'`
    outfile="${outdir}/${prefix}.all_covariates.txt"
    Rscript cov-format.R ${outdir}/${file} ${outfile} &
done

wait

# ---------------------------------------------------------------------------

base="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8"
outdir="/mnt/lab_data/montgomery/nicolerg/eqtl"

genotype_pcs="${outdir}_supporting/genot_5pcs.txt" # header must be "ID" "GTEX###"...

source activate py35

for prefix in `ls ${outdir} | grep "all_covariates" | sed 's/\..*//'`
do
	taskset -c 12-47 python /users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/combine_covariates.py ${outdir}/${prefix}.all_covariates.txt ${prefix} --genotype_pcs ${genotype_pcs} -o ${outdir}
	# This step generates a combined covariates file, 
	# containing genotype PCs, PEER factors, and additional explicit covariates (e.g., genotyping platform)
	# The covariate files should have one covariate per row, with an identifier in the first column, and a header line with sample identifiers
	# This step will generate the file ${prefix}.combined_covariates.txt
done

source deactivate

# ---------------------------------------------------------------------------

source activate py35

# for prefix in `ls ${outdir} | grep "combined_covariates" | sed 's/\..*//'`
# do
# 	python /users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} --covariates ${prefix}.combined_covariates.txt --permute 1000 10000 --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 30 --threads 36 -o ${outdir} &
# 	# The following files will be generated:
# 	# ${prefix}.allpairs.txt.gz
# 	# ${prefix}.egenes.txt.gz
# done

prefix="Whole_Blood"
python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py ${vcf} ${outdir}/${prefix}.expression.bed.gz ${prefix} --covariates ${outdir}/${prefix}.combined_covariates.txt --permute 1000 10000 --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 100 --threads 36 -o ${outdir}

source deactivate