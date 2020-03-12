tissue="Whole_Blood"
chr="1"
outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood-new"
indir="${outdir}/chunks"

geno="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/gtex.MAF01.phased.chr${chr}.vcf.gz"
localcov="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr${chr}.hg19.localcov.tsv.gz"
globalcov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/${tissue}.v8.covariates.txt"
pheno="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

# split expression into chunks
exprdir="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices"
exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz

exp="${outdir}/${tissue}.chr${chr}.expression.bed"
if [ ! -d "$indir" ]; then
	mkdir ${indir}
fi
zgrep -E "^chr${chr}	" ${exprfile} > ${exp}
split -d -l 30 ${exp} "${indir}/${tissue}.chr${chr}.expression."

header=${indir}/header.txt
zcat ${exprfile} | head -1 > ${header}
for file in `ls ${indir} | grep "expression"`
do
	cat ${header} > tmp.txt
	cat ${indir}/${file} >> tmp.txt
	rm ${indir}/${file}
	mv tmp.txt ${indir}/${file}
done
gzip ${indir}/${tissue}.chr${chr}.expression*

for chunk in `ls ${indir} | grep "chr${chr}.expression"`
do
	gene=${indir}/${chunk}
	suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`
	echo "Starting chr $chr chunk $suf"
	out="${outdir}/${tissue}-LAVA-global-allpairs-chr${chr}-${suf}.tsv"
	taskset -c 43-47 Rscript ~/gtex-admix/scripts/geuvadis/lava-lcl-afr.R ${geno} ${gene} ${globalcov} ${localcov} ${pheno} ${out}
	
	exit

	running=`ps -ef | grep "nicolerg" | grep "lava-lcl-afr" | wc -l`
	while [ $running -gt 9 ] # run max 10 jobs at a time
	do
		sleep 180
		running=`ps -ef | grep "nicolerg" | grep "lava-lcl-afr" | wc -l`
	done
done
