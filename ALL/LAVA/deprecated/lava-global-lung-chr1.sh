tissue="Lung"
base="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed"

exprdir="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices"
exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz
globalcov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/${tissue}.v8.covariates.txt"
pheno="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
admix_ids="/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt"

# set up directories
indir=${base}/${tissue}
if [ ! -d "$indir" ]; then
	mkdir $indir
done
edir=${indir}/chunks
if [ ! -d "$edir" ]; then
	mkdir $edir
done
for chr in {1..22}; do
	outdir=${indir}/chr${chr}
	if [ ! -d "$outdir" ]; then
		mkdir $outdir
	done
done
vcfdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno"
if [ ! -d "$vcfdir" ]; then
	mkdir $vcfdir
done

# filter VCFs by AC (only do once per group of subjects)

vcfin="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19"

if [ ! -e "$vcfdir/gtex.admixed.MAC10.phased.chr${chr}.recode.vcf.gz" ]; then
	for chr in {1..22}; do
		vcf_unfilt=${vcfin}/gtex.MAF01.phased.chr${chr}.vcf.gz
		vcf_out=${vcfdir}/gtex.admixed.MAC10.phased.chr${chr}
		taskset -c 24-47 vcftools \
			--gzvcf ${vcf_unfilt} \
			--keep ${admix_ids} \
			--recode-INFO-all \
			--mac 10 \
			--max-alleles 2 \
			--min-alleles 2 \
			--out ${vcf_out} \
			--recode &
	done

	wait

	gzip ${vcfdir}/gtex.admixed.MAC10.phased.chr*

fi 

for chr in {1..22}; do

	geno="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/gtex.MAF01.phased.chr${chr}.vcf.gz"
	localcov="/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr${chr}.hg19.localcov.tsv.gz"

	outdir="${indir}/chr${chr}"

	# split expression into chunks
	exp="${outdir}/${tissue}.chr${chr}.expression.bed"
	zgrep -E "^chr${chr}	" ${exprfile} > ${exp}
	split -d -l 50 ${exp} "${edir}/${tissue}.chr${chr}.expression."

	header=${edir}/header.chr${chr}.txt
	zcat ${exprfile} | head -1 > ${header}
	for file in `ls ${edir} | grep "chr${chr}.expression"`
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
		out="${outdir}/${tissue}-LAVA-global-allpairs-chr${chr}-${suf}.tsv"
		taskset -c 40-47 Rscript ~/gtex-admix/scripts/LAVA/lava-global-${tissue}-chr1.R ${chr} ${tissue} ${gene} ${out} &

		running=`ps -ef | grep "nicolerg" | grep "lava-global-${tissue}" | wc -l`
		while [ $running -gt 8 ] # run max 10 jobs at a time
		do
			sleep 180
			running=`ps -ef | grep "nicolerg" | grep "lava-global-${tissue}" | wc -l`
		done
		
	done

	wait

	# merge chunks 
	f=`ls ${outdir} | grep "${tissue}-LAVA-global-allpairs" | head -1`
	head -1 ${outdir}/${f} > ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv
	for chunk in `ls ${outdir} | grep "${tissue}-LAVA-global-allpairs-chr${chr}"`; do
		file=${outdir}/${chunk}
		cat ${file} | sed -e '1d' >> ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv
		rm ${file}
	done
	gzip ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv

done

# merge into whole genome  
zcat ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr1.tsv.gz | head -1 > ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv
for file in ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr{1..22}.tsv.gz; do
	zcat $file | sed -e '1d' >> ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv
done
gzip ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv

# # define eQTL genes
# python /users/nicolerg/gtex-admix/scripts/LAVA/extract-egene-lava.py ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv.gz
# Rscript /users/nicolerg/gtex-admix/scripts/LAVA/pval-correct-lava.R ${outdir}

# # assuming GEU gene sets have already been defined, parse regulatory enrichment
# geudir="/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl"
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py ${geudir}/eur.eqtl.genes.filt.txt 'eur.eqtl.filt'
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py ${geudir}/afr.eqtl.genes.txt 'afr.eqtl'
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py ${geudir}/eur.eqtl.genes.txt 'eur.eqtl'
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py

# # plot regulatory enrichment
# regdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/${tissue}/aa-only/reg-enrich"
# Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${regdir}/global.regenrich.afr.eqtl.tsv.gz ${regdir}/lava.regenrich.afr.eqtl.tsv.gz 'afr.eqtl'
# Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${regdir}/global.regenrich.eur.eqtl.tsv.gz ${regdir}/lava.regenrich.eur.eqtl.tsv.gz 'eur.eqtl'
# Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${regdir}/global.regenrich.eur.eqtl.filt.tsv.gz ${regdir}/lava.regenrich.eur.eqtl.filt.tsv.gz 'eur.eqtl.filt'
# Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${regdir}/global.regenrich.all.eqtl.genes.tsv.gz ${regdir}/lava.regenrich.all.eqtl.genes.tsv.gz 'all.eqtl.genes'
