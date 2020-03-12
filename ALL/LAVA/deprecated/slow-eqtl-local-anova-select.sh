globalcov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Lung.v8.covariates.txt"
pheno="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/anova"
indir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/anova/filtered/chunks"

geno="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/anova/filtered/chr1-geno-filt-byafrsig.tsv"
localcov="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/anova/filtered/chr1-localcov-filt-byafrsig.tsv"

for chunk in `ls $indir | grep "chr1.expression"`
do
	gene=${indir}/${chunk}
	suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`
	echo "Starting chunk $suf"
	out1="${outdir}/lung-local-allpairs-chr1-${suf}.tsv"
	taskset -c 43-47 Rscript /users/nicolerg/gtex-admix/scripts/local_eqtl/slow-eqtl-local-anova-select.R ${geno} ${gene} ${globalcov} ${localcov} ${pheno} ${out1} &
	
	running=`ps -ef | grep "nicolerg" | grep "slow-eqtl-local" | wc -l`
	while [ $running -gt 8 ] # run max 8 jobs at a time
	do
		sleep 180
		running=`ps -ef | grep "nicolerg" | grep "slow-eqtl-local" | wc -l`
	done
done
