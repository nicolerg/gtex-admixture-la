globalcov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Lung.v8.covariates.txt"
pheno="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/anova"
indir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chunks"

for chr in 1
do	
	geno="/mnt/lab_data/montgomery/nicolerg/filtered/gtex.MAF01.phased.chr${chr}.vcf.gz"
	localcov="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr${chr}.hg19.localcov.tsv.gz"
	for chunk in `ls $indir | grep "chr${chr}.expression"`
	#for chunk in ${indir}/lung.chr1.expression.90{24..71}.gz
	do
		gene=${indir}/${chunk}
		suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`
		echo "Starting chr $chr chunk $suf"
		out1="${outdir}/lung-local-allpairs-chr${chr}-${suf}.tsv"
		taskset -c 43-47 Rscript slow-eqtl-local-optimized-comparison.R ${geno} ${gene} ${globalcov} ${localcov} ${pheno} ${out1} &
		
		running=`ps -ef | grep "nicolerg" | grep "slow-eqtl-local" | wc -l`
		while [ $running -gt 5 ] # run max 5 jobs at a time
		do
			sleep 180
			running=`ps -ef | grep "nicolerg" | grep "slow-eqtl-local" | wc -l`
		done
	done
done

## for a single chromosome, wholeblood:

# indir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets"

# geno="/mnt/lab_data/montgomery/nicolerg/filtered/gtex.phased.MAF01.hg19.chr21.vcf"
# globalcov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Adipose_Subcutaneous.v8.covariates.txt"
# localcov="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/chr21.hg19.localcov.tsv"
# pheno="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
# outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood"

# for chunk in {1..11}
# do 
# 	gene="${indir}/chunk${chunk}.bed"
# 	#outfile="${outdir}/wholeblood-global-chr21-chunk${chunk}.tsv"
# 	#taskset -c 0-24 Rscript slow-eqtl-global.R ${geno} ${gene} ${globalcov} ${localcov} ${pheno} ${outfile} &
# 	out1="${outdir}/wholeblood-local-chr21-chunk${chunk}.tsv"
# 	out2="${outdir}/wholeblood-local-unkcnt-chr21-chunk${chunk}.tsv"
# 	taskset -c 0-24 Rscript slow-eqtl-local.R ${geno} ${gene} ${globalcov} ${localcov} ${pheno} ${out1} ${out2} &
# done