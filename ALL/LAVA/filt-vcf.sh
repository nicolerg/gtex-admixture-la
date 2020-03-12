# filter VCFs for custom GTEx eQTL calling

indir='/mnt/lab_data/montgomery/nicolerg/genotypes-hg19'
outdir=${indir}/gtex-admixed-hg19

# for chr in {1..22}; do
# 	infile=${indir}/gtex.MAF01.phased.chr${chr}.vcf.gz
# 	outprefix=${outdir}/gtex.admixed.maf.chr${chr}
# 	subjs=${indir}/lung-admix-subj.txt

# 	vcftools --gzvcf ${infile} \
# 		--out ${outprefix} \
# 		--keep ${subjs} \
# 		--maf 0.02 \
# 		--min-alleles 2 \
# 		--max-alleles 2 \
# 		--recode \
# 		--recode-INFO-all & # keep all original INFO in output file
# done

# wait

# for chr in {1..22}; do
# 	bgzip ${outdir}/gtex.admixed.maf.chr${chr}.recode.vcf &
# done

# wait

# now filter chromosome-level allpairs files by filtered VCFs
for chr in {1..22}; do 
	python ~/gtex-admix/scripts/LAVA/filt-allpairs-snps.py ${chr} &

	running=`ps -ef | grep "nicolerg" | grep "filt-allpairs-snps" | wc -l`
	while [ $running -gt 8 ] # run max 10 jobs at a time
	do
		sleep 60
		running=`ps -ef | grep "nicolerg" | grep "filt-allpairs-snps" | wc -l`
	done
done
