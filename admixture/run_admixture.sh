# run ADMIXTURE on VCFs 
# following http://software.genetics.ucla.edu/admixture/admixture-manual.pdf

set -e

indir=/projects/nicolerg/gtex-tg-merged
outdir=/projects/nicolerg/gtex-tg-merged/admixture_20190819-3
tg_meta=1kg-subject-popcodes.tsv
admix=gtex-admixed0.9.txt

# ${indir} contains chromosome-level VCFs generated as follows, using Snakefile:
# 1. GTEx and 1000 Genomes filtered with MAF > 0.05
# 2. find intersection of each chromosome
# 3. merge intersection

# now do the following for ADMIXTURE:
# 1. filter subject (117 admixed individuals and 1000 Genomes reference panels) and apply MAF filter of 0.10
# 2. remove SNPs with missing values
# 3. LD-prune each chromosome
# 4. concatenate chromosomes
# 5. generate .pop file 
# 6. run ADMIXTURE 

# run admixture on 117 GTEx individuals with any set of reference panels

mkdir -p ${outdir}
cd ${outdir}

# select reference populations 
vcf_header=${indir}/genotyped_subjects.txt
zgrep -m1 -E "^#CHROM" ${indir}/merged_hg38_chr22.vcf.gz | cut --complement -f1-9 | sed "s/	/\n/g" > ${vcf_header}

# string with user-defined SUPERPOPULATIONS
# modify this Rscript to adjust how 1000 Genomes reference panels are selected
Rscript /scripts/gtex-admixture-la/admixture/select_reference_populations.R 'EUR,EAS,AFR' ${vcf_header} ${tg_meta} ${admix}
K=3 # number of superpopulations 

cut -f1 keep_subjects.txt > bcftools_keep.txt

format_chrom () {
	local chr=$1
	# filter subjects and re-apply MAF filter
	bcftools view \
		--force-samples \
		-S bcftools_keep.txt \
		-Oz \
		--min-af 0.1 \
		-o merged_hg38_chr${chr}_filt_MAF10.vcf.gz \
		../merged_hg38_chr${chr}.vcf.gz 
	
	# remove SNPs with any missing values 
	bcftools view \
		-g ^miss \
		-Oz \
		-o merged_hg38_chr${chr}_filt_MAF10_nomissing.vcf.gz \
		merged_hg38_chr${chr}_filt_MAF10.vcf.gz 
	
	#index
	tabix -p vcf merged_hg38_chr${chr}_filt_MAF10_nomissing.vcf.gz 
	
	# convert to binary .BED
	plink --vcf merged_hg38_chr${chr}_filt_MAF10_nomissing.vcf.gz \
		--make-bed \
		--vcf-half-call m \
		--out merged_hg38_chr${chr}_filt_MAF10_nomissing 
	
	# LD-prune each chromosome 
	plink --bfile merged_hg38_chr${chr}_filt_MAF10_nomissing \
		--indep-pairwise 40 10 0.1 \
		--out merged_hg38_chr${chr}_filt_MAF10_nomissing 
	plink --bfile merged_hg38_chr${chr}_filt_MAF10_nomissing \
		--extract merged_hg38_chr${chr}_filt_MAF10_nomissing.prune.in \
		--make-bed \
		--out merged_hg38_chr${chr}_filt_MAF10_nomissing_pruned 
}

# this will occupy 22 cores 
for chr in {1..22}; do format_chrom "${chr}" & done
wait

# concatenate all chromosomes 
ls | grep "pruned.bed" > bed
ls | grep "pruned.bim" > bim
ls | grep "pruned.fam" > fam
paste bed bim fam > files_to_merge.txt
sed -i -e '1d' files_to_merge.txt

plink --bfile merged_hg38_chr${chr}_filt_MAF10_nomissing_pruned \
	--merge-list files_to_merge.txt \
	--make-bed \
	--out merged_hg38_filt_MAF10_nomissing_pruned

# remove intermediate files 
rm merged_hg38_chr*_filt_MAF10_nomissing.* \
	merged_hg38_chr*_filt_MAF10_nomissing_pruned.* \
	bed bim fam 

# generate .pop file 
# I want to keep all individuals in gtex-admixed0.9.txt
# as well as certain 1000 Genomes populations. 
# we need to do multiple iterations with different reference panels, but include the following to start:
# Asian (N = 207)
#   CHB (Han Chinese in Beijing, China)
#   JPT (Japanese in Tokyo, Japan)
# African (N = 405)
#   YRI (Yoruba in Ibadan, Nigeria)
#   GWD (Gambian in Western Divisions in the Gambia)
#   MSL (Mende in Sierra Leone)
#   ESN (Esan in Nigeria)
# European (N = 99)
#   CEU
# (same as RFMix reference panels)
tg=/projects/1KG/1kg-subject-popcodes.tsv

Rscript /scripts/gtex-admixture-la/admixture/order_pop_file.R merged_hg38_filt_MAF10_nomissing_pruned.fam ${tg} 
cut -f3 merged_hg38_filt_MAF10_nomissing_pruned.pop.txt > merged_hg38_filt_MAF10_nomissing_pruned.pop

# add software to path 
export PATH=$PATH:/users/nicolerg/admixture_linux-1.3.0
# run admixture 
#K=3 # number of reference populations included in the .pop file 
admixture -j12 merged_hg38_filt_MAF10_nomissing_pruned.bed ${K} --supervised

paste *pop.txt *.Q | sed "s/ /	/g" > merged_hg38_filt_MAF10_nomissing_pruned.${K}.Q.subj

