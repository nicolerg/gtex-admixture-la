tissue="Lung"
outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/fastqtl"
gtexdir="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8"
globalcov=${gtexdir}/eqtl/GTEx_Analysis_v8_eQTL_covariates/${tissue}.v8.covariates.txt
pheno=${gtexdir}/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt

exprdir=${gtexdir}/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices
exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz

exprout="$outdir/expr"
if [ ! -d "$outdir/expr" ]; then
	mkdir $outdir/expr
fi
if [ ! -d "$outdir/geno" ]; then
	mkdir $outdir/geno
fi
if [ ! -d "$outdir/cov" ]; then
	mkdir $outdir/cov
fi
if [ ! -d "$outdir/eqtl" ]; then
	mkdir $outdir/eqtl
fi

# make covariates files (30 PEER)
l=`grep -n "InferredCov30" ${globalcov} | sed "s/:.*//"`
head -$l ${globalcov} > $outdir/cov/${tissue}.aa.cov.txt
cov=$outdir/cov/${tissue}.aa.cov.txt

# get list of admixed individuals 
cut -f1,5 ${pheno} | grep -E "	2|	1" | cut -f1 > ${outdir}/geno/aa-subjs.txt

if [ ! -e "${exprout}/gtex.lung.admix.expression.bed" ]; then
	# filter expression by admixed 
	zcat ${exprfile} | cut -f1,2,3,4 > ${exprout}/gtex.lung.admix.expression.bed
	for sub in `cat ${outdir}/geno/aa-subjs.txt`; do
		col=`zcat ${exprfile} | head -1 | sed 's/	/\n/g' | grep -n "$sub" | sed "s/:.*//"`
		if [ "$col" != "" ]; then
			zcat ${exprfile} | cut -f$col > ${exprout}/tmp
			paste ${exprout}/gtex.lung.admix.expression.bed ${exprout}/tmp > ${exprout}/tmp.2
			rm ${exprout}/gtex.lung.admix.expression.bed
			mv ${exprout}/tmp.2 ${exprout}/gtex.lung.admix.expression.bed
		fi
	done
fi	
filtexpr=${exprout}/gtex.lung.admix.expression.bed

# get list of admixed individuals with expression
cut --complement -f1,2,3,4 ${filtexpr} | head -1 | sed 's/	/\n/g' > ${outdir}/geno/aa-filt.txt

# filter VCFs by subject (admixed with expression) and MAF
vcfdir="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/gtex-phased-MAF0-hg19"
if [ ! -e "${outdir}/geno/gtex.admix.chr22.recode.vcf.gz" ]; then
	for chr in {1..22}; do 
		vcf=${vcfdir}/gtex.chr${chr}.vcf.gz
		vcftools --gzvcf ${vcf} --out $outdir/geno/gtex.admix.chr${chr} --keep ${outdir}/geno/aa-filt.txt --recode --recode-INFO-all --maf 0.05 --min-alleles 2 --max-alleles 2 &
	done
	wait
	for chr in {1..22}; do 
		bgzip -f $outdir/geno/gtex.admix.chr${chr}.recode.vcf
		tabix -f -p vcf $outdir/geno/gtex.admix.chr${chr}.recode.vcf.gz
	done
fi

# run FastQTL
for chr in {1..22}; do

	vcf=$outdir/geno/gtex.admix.chr${chr}.recode.vcf.gz
	nchunks=50
	
	# select only expression in this chromosome
	head -1 "${filtexpr}" > "${exprout}/gtex.admix.chr${chr}.expr.bed"
	grep -E "^chr${chr}	" ${filtexpr} >> "${exprout}/gtex.admix.chr${chr}.expr.bed"
	bgzip "${exprout}/gtex.admix.chr${chr}.expr.bed"
	tabix -p bed "${exprout}/gtex.admix.chr${chr}.expr.bed.gz"
	expr="${exprout}/gtex.admix.chr${chr}.expr.bed.gz"

	# permutation pass
	for j in $(seq 1 $nchunks); do
		taskset -c 30-47 /users/nicolerg/FastQTL/bin/fastQTL \
		--vcf ${vcf} \
		--bed ${expr} \
		--window 1e6 \
		--out ${outdir}/eqtl/gtex.admix.permuted.chr${chr}.chunk${j}.txt.gz \
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

# remove empty files
gunzip *.gz
rm *.gz 

# merge each chr
for chr in {1..22}; do 
	echo "gene_id num_var shape1_mle shape2_mle dummy variant tss_distance pval_nominal slope pval_perm pval_beta" > ${outdir}/eqtl/gtex.admix.permuted.chr${chr}.txt
	for j in $(seq 1 $nchunks); do
		file=${outdir}/eqtl/gtex.admix.permuted.chr${chr}.chunk${j}.txt
		if [ -e "$file" ]; then 
			cat $file >> ${outdir}/eqtl/gtex.admix.permuted.chr${chr}.txt
		fi
		# rm $file 
	done
	gzip ${outdir}/eqtl/gtex.admix.permuted.chr${chr}.txt
done

# merge chr into whole genome
if [ ! -e "${outdir}/eqtl/gtex.admix.permuted.txt.gz" ]; then
	zcat ${outdir}/eqtl/gtex.admix.permuted.chr1.txt.gz | head -1 > ${outdir}/eqtl/gtex.admix.permuted.txt
	for file in ${outdir}/eqtl/gtex.admix.permuted.chr{1..22}.txt.gz; do 
		zcat $file | sed -e '1d' >> ${outdir}/eqtl/gtex.admix.permuted.txt
		# rm $file
	done 
	gzip ${outdir}/eqtl/gtex.admix.permuted.txt
fi
