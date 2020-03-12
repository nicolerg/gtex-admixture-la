chr=$1
pop=$2

if [ "$pop" -eq "EUR" ]; then 

	# for EUR into AFR
	indir="/mnt/lab_data/montgomery/nicolerg/geuvadis"
	vcfdir=${indir}/genotype/EUR-into-AFR
	vcfin="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/1000genomes-phased-hg19"

	vcf=${vcfin}/tg.chr${chr}.vcf.gz

	eur_vcf=${vcfdir}/tg.EUR.chr${chr}.recode.vcf.gz
	afr_vcf=${vcfdir}/tg.AFR.chr${chr}.recode.vcf.gz

	zcat "${indir}/expression/GV-AFR.expression.bed.gz" | cut --complement -f1,2,3,4 | head -1 | sed 's/	/\n/g' > ${vcfdir}/afr-subjs.txt
	zcat "${indir}/expression/GV-EUR.expression.bed.gz" | cut --complement -f1,2,3,4 | head -1 | sed 's/	/\n/g' > ${vcfdir}/eur-subjs.txt

	if [ ! -e "$afr_vcf" -o ! -e "$eur_vcf" ]; then

		# EUR
		vcftools --gzvcf ${vcf} --out ${vcfdir}/tg.EUR.chr${chr} --keep ${vcfdir}/eur-subjs.txt --recode --recode-INFO-all --maf 0.05 --min-alleles 2 --max-alleles 2
		bgzip -f ${vcfdir}/tg.EUR.chr${chr}.recode.vcf
		tabix -f -p vcf ${vcfdir}/tg.EUR.chr${chr}.recode.vcf.gz
		zcat ${vcfdir}/tg.EUR.chr${chr}.recode.vcf.gz | cut -f3 | sed -e '1d' > ${vcfdir}/keep-snp-chr${chr}
		#AFR 
		vcftools --gzvcf ${vcf} --out ${vcfdir}/tg.AFR.chr${chr} --keep ${vcfdir}/afr-subjs.txt --recode --recode-INFO-all --snps ${vcfdir}/keep-snp-chr${chr} --min-alleles 2 --max-alleles 2
		bgzip -f ${vcfdir}/tg.AFR.chr${chr}.recode.vcf
		tabix -f -p vcf ${vcfdir}/tg.AFR.chr${chr}.recode.vcf.gz
	fi
fi

if [ "$pop" -eq "AFR" ]; then 

	# for AFR into EUR 
	indir="/mnt/lab_data/montgomery/nicolerg/geuvadis"
	vcfdir=${indir}/genotype/AFR-into-EUR
	vcfin="/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/1000genomes-phased-hg19"

	vcf=${vcfin}/tg.chr${chr}.vcf.gz

	eur_vcf=${vcfdir}/tg.EUR.chr${chr}.recode.vcf.gz
	afr_vcf=${vcfdir}/tg.AFR.chr${chr}.recode.vcf.gz

	zcat "${indir}/expression/GV-AFR.expression.bed.gz" | cut --complement -f1,2,3,4 | head -1 | sed 's/	/\n/g' > ${vcfdir}/afr-subjs.txt
	zcat "${indir}/expression/GV-EUR.expression.bed.gz" | cut --complement -f1,2,3,4 | head -1 | sed 's/	/\n/g' > ${vcfdir}/eur-subjs.txt

	# AFR
	vcftools --gzvcf ${vcf} --out ${vcfdir}/tg.AFR.chr${chr} --keep ${vcfdir}/afr-subjs.txt --recode --recode-INFO-all --maf 0.05 --min-alleles 2 --max-alleles 2
	bgzip -f ${vcfdir}/tg.AFR.chr${chr}.recode.vcf
	tabix -f -p vcf ${vcfdir}/tg.AFR.chr${chr}.recode.vcf.gz
	zcat ${vcfdir}/tg.AFR.chr${chr}.recode.vcf.gz | cut -f3 | sed -e '1d' > ${vcfdir}/keep-snp-chr${chr}
	# EUR 
	vcftools --gzvcf ${vcf} --out ${vcfdir}/tg.EUR.chr${chr} --keep ${vcfdir}/eur-subjs.txt --recode --recode-INFO-all --snps ${vcfdir}/keep-snp-chr${chr} --min-alleles 2 --max-alleles 2
	bgzip -f ${vcfdir}/tg.EUR.chr${chr}.recode.vcf
	tabix -f -p vcf ${vcfdir}/tg.EUR.chr${chr}.recode.vcf.gz
fi
