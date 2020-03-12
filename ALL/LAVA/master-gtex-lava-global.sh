# DEPRECATED: see `batch-master-lava-global.sh`


# VCF filtering is the only important part of this script:

# filter VCFs by AC (only do once per group of subjects)

vcfin="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19"
vcfdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno"

if [ ! -e "$vcfdir/gtex.admixed.MAC10.phased.chr1.recode.vcf.gz" ]; then
	for chr in {1..22}; do
		vcf_unfilt=${vcfin}/gtex.MAF01.phased.chr${chr}.vcf.gz
		vcf_out=${vcfdir}/gtex.admixed.MAC10.phased.chr${chr}
		taskset -c 24-47 vcftools \
			--gzvcf ${vcf_unfilt} \
			--keep ${admix_ids} \
			--recode-INFO-all \
			--maf 0.01 \
			--mac 10 \
			--max-alleles 2 \
			--min-alleles 2 \
			--out ${vcf_out} \
			--recode &
	done

	wait

	gzip ${vcfdir}/gtex.admixed.MAC10.phased.chr*

fi 

# for chr in {1..22}; do

# 	for chunk in `ls ${edir} | grep "chr${chr}.expression"`; do
# 		gene=${edir}/${chunk}
# 		suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`
# 		echo "Starting chr $chr chunk $suf"
# 		out="${outdir}/${tissue}-LAVA-global-allpairs-chr${chr}-${suf}.tsv"
# 		taskset -c 42-47 Rscript ${srcdir}/gtex-lava-global.R ${chr} ${tissue} ${gene} ${globalcov} ${out} &

# 		running=`ps -ef | grep "nicolerg" | grep "gtex-lava-global.R" | wc -l`
# 		while [ $running -gt 5 ]; do
# 			sleep 180
# 			running=`ps -ef | grep "nicolerg" | grep "gtex-lava-global.R" | wc -l`
# 		done
		
# 	done

# 	wait

# 	# merge chunks 
# 	f=`ls ${outdir} | grep "${tissue}-LAVA-global-allpairs" | head -1`
# 	head -1 ${outdir}/${f} > ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv
# 	for chunk in `ls ${outdir} | grep "${tissue}-LAVA-global-allpairs-chr${chr}"`; do
# 		file=${outdir}/${chunk}
# 		cat ${file} | sed -e '1d' >> ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv
# 		rm ${file}
# 	done
# 	gzip ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv

# done

# # merge into whole genome  
# outdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}

# zcat ${outdir}/chr1/${tissue}-LAVA-global-merged-allpairs-chr1.tsv.gz | head -1 > ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv
# for dir in chr{1..22}; do
# 	file=${outdir}/${dir}/${tissue}-LAVA-global-merged-allpairs-${dir}.tsv.gz
# 	zcat $file | sed -e '1d' >> ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv
# done
# gzip ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv

# =============================

# # summarize eGenes; define eQTL genes
# python ${srcdir}/extract-egene-lava.py ${indir}/${tissue}-LAVA-global-merged-allpairs.tsv.gz

# # parse regulatory enrichment for different gene sets
# diff_egenes='/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/top-100-diff-by-effect.txt'
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py ${diff_egenes} 'top100.diff'
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py



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
