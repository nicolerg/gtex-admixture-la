tissue="Lung"
indir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chunks"
globalcov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/${tissue}.v8.covariates.txt"
pheno="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

exprdir="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices"
exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz

chr='1'

geno="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/gtex.MAF01.phased.chr${chr}.vcf.gz"
localcov="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr${chr}.hg19.localcov.tsv.gz"

outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/chr${chr}"
if [ ! -d "$outdir" ]; then
	mkdir ${outdir}
fi

# split expression into chunks
exp="${outdir}/${tissue}.chr${chr}.expression.bed"
if [ ! -d "$indir" ]; then
	mkdir ${indir}
fi
zgrep -E "^chr${chr}	" ${exprfile} > ${exp}
split -d -l 50 ${exp} "${indir}/${tissue}.chr${chr}.expression."

header=${indir}/header.chr${chr}.txt
zcat ${exprfile} | head -1 > ${header}
for file in `ls ${indir} | grep "chr${chr}.expression"`
do
	cat ${header} > tmp.chr${chr}.txt
	cat ${indir}/${file} >> tmp.chr${chr}.txt
	rm ${indir}/${file}
	mv tmp.chr${chr}.txt ${indir}/${file}
done
gzip ${indir}/${tissue}.chr${chr}.expression*

for chunk in `ls ${indir} | grep "chr${chr}.expression"`; do
	gene=${indir}/${chunk}
	suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`
	echo "Starting chr $chr chunk $suf"
	out="${outdir}/lung-LAVA-global-allpairs-chr${chr}-${suf}.tsv"
	taskset -c 40-47 Rscript ~/gtex-admix/scripts/LAVA/lava-global-lung-test.R ${chr} ${tissue} ${gene} ${out} 

	exit
done
