# plot expression per ${poly}opean genotype for ${poly} eQTLs that have 0 effect size in YRI (GEU eQTL calling)

# eqtldir="/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR"
# vcfdir="/mnt/lab_data/montgomery/nicolerg/geuvadis/genotype/AFR-into-EUR"
# mono="EUR"
# poly="AFR"

# if [ ! -e "${eqtldir}/${poly}.eqtl.${mono}.0slope.txt" ]; then 
# 	Rscript ~/gtex-admix/scripts/geuvadis/geu-custom-effect-sizes.R ${eqtldir}
# fi

# if [ ! -e "${vcfdir}/tg.${mono}.recode.vcf.gz" ]; then 

# 	zcat ${vcfdir}/tg.${mono}.chr1.recode.vcf.gz | head -1 > ${vcfdir}/tg.${mono}.recode.vcf
# 	for file in ${vcfdir}/tg.${mono}.chr{1..22}.recode.vcf.gz; do
# 		zcat $file | sed -e '1d' >> ${vcfdir}/tg.${mono}.recode.vcf
# 	done
# 	bgzip ${vcfdir}/tg.${mono}.recode.vcf
# fi

# if [ ! -e "${vcfdir}/tg.${poly}.recode.vcf.gz" ]; then 

# 	zcat ${vcfdir}/tg.${poly}.chr1.recode.vcf.gz | head -1 > ${vcfdir}/tg.${poly}.recode.vcf
# 	for file in ${vcfdir}/tg.${poly}.chr{1..22}.recode.vcf.gz; do
# 		zcat $file | sed -e '1d' >> ${vcfdir}/tg.${poly}.recode.vcf
# 	done
# 	bgzip ${vcfdir}/tg.${poly}.recode.vcf
# fi

# cut -f2 ${eqtldir}/${poly}.eqtl.${mono}.0slope.txt | sed -e '1d' > ${eqtldir}/${mono}.0slope.${poly}.eqtl.snps

# # filter AFR vcf by these SNPs

# if [ ! -e "${vcfdir}/pop-diff/tg.${mono}.recode.0slope.recode.vcf.gz" ]; then 
# 	vcftools --gzvcf ${vcfdir}/tg.${mono}.recode.vcf.gz \
# 		--recode-INFO-all \
# 		--out ${vcfdir}/pop-diff/tg.${mono}.recode.0slope \
# 		--recode \
# 		--snps ${eqtldir}/${mono}.0slope.${poly}.eqtl.snps

# 		# get freq for these SNPs

# 		vcftools --vcf ${vcfdir}/pop-diff/tg.${mono}.recode.0slope.recode.vcf \
# 			--freq \
# 			--out ${vcfdir}/pop-diff/tg.${mono}.0slope

# fi

# zcat ${vcfdir}/pop-diff/tg.${mono}.recode.0slope.recode.vcf.gz | cut -f3 > ${vcfdir}/pop-diff/rsid
# paste ${vcfdir}/pop-diff/tg.${mono}.0slope.frq ${vcfdir}/pop-diff/rsid > ${vcfdir}/pop-diff/tg.${mono}.0slope.freq
# bgzip ${vcfdir}/pop-diff/tg.${mono}.recode.0slope.recode.vcf

# if [ ! -e "${vcfdir}/pop-diff/tg.${poly}.recode.0slope.recode.vcf.gz" ]; then 

# 	vcftools --gzvcf ${vcfdir}/tg.${poly}.recode.vcf.gz \
# 		--recode-INFO-all \
# 		--out ${vcfdir}/pop-diff/tg.${poly}.recode.0slope \
# 		--recode \
# 		--snps ${eqtldir}/${mono}.0slope.${poly}.eqtl.snps
# 	bgzip ${vcfdir}/pop-diff/tg.${poly}.recode.0slope.recode.vcf

# 	# vcftools --gzvcf ${vcfdir}/pop-diff/tg.${poly}.recode.0slope.recode.vcf.gz \
# 	# 	--freq \
# 	# 	--out tg.${poly}.0slope
# fi

# python ~/gtex-admix/scripts/geuvadis/eqtl-snp-expression.py ${eqtldir}/${poly}.eqtl.${mono}.0slope.txt ${vcfdir}/pop-diff/tg.${mono}.0slope.freq ${vcfdir}/pop-diff/tg.${poly}.recode.0slope.recode.vcf.gz ${mono}
# Rscript ~/gtex-admix/scripts/geuvadis/plot-expr-by-genotype.R ${eqtldir}/${poly}.eqtl.0${mono}.var.expr.tsv.gz ${eqtldir}
# python ~/gtex-admix/scripts/geuvadis/parse-tstat.py ${eqtldir}/GEU.${poly}.eqtl.0${mono}.tstat.tsv ${eqtldir}
# Rscript ~/gtex-admix/scripts/geuvadis/tstat-boxplot.R 

# now I want to look at expression for SNPs that are eQTLs in both populations ===========================================
vcfbase="/mnt/lab_data/montgomery/nicolerg/geuvadis/genotype"
base="/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl"

# sed -e '1d' ${base}/EUR-into-AFR/EUR.geu.egenes.5fdr.txt | cut -f1,6 > ${base}/eur #EUR eQTLs
# sed -e '1d' ${base}/AFR-into-EUR/AFR.geu.egenes.5fdr.txt | cut -f1,6 > ${base}/afr #AFR eQTLs
# sort ${base}/afr ${base}/eur | uniq -d > ${base}/overlapping/overlapping.eqtls
# cut -f2 ${base}/overlapping/overlapping.eqtls > ${base}/overlapping/overlapping.snps

# # filter AFR vcf by these SNPs

# vcftools --gzvcf ${vcfbase}/AFR-into-EUR/tg.AFR.recode.vcf.gz \
# 	--recode-INFO-all \
# 	--out ${base}/overlapping/tg.AFR.recode.overlap \
# 	--recode \
# 	--snps ${base}/overlapping.snps &

# # filter EUR vcf by these SNPs

# vcftools --gzvcf ${vcfbase}/EUR-into-AFR/tg.EUR.recode.vcf.gz \
# 	--recode-INFO-all \
# 	--out ${base}/overlapping/tg.EUR.recode.overlap \
# 	--recode \
# 	--snps ${base}/overlapping.snps &

# wait

# bgzip ${base}/overlapping/tg.EUR.recode.overlap.recode.vcf
# bgzip ${base}/overlapping/tg.AFR.recode.overlap.recode.vcf

python ~/gtex-admix/scripts/geuvadis/overlap-eqtl-snp-expression.py
Rscript ~/gtex-admix/scripts/geuvadis/overlap-plot-expr-by-genotype.R 
