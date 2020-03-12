#!/bin/bash

set -e

################################################################################################################
## Generate merged VCFs (GTEx and 1000 Genomes, by chromosome)
################################################################################################################

outdir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged

# 1000 Genomes VCFs, already split by chromosome:
tg=/mnt/lab_data/montgomery/shared/1KG/hg38

# GTEx VCF, needs to be split by chromosome:
gtex_vcf=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
# split GTEx VCF by chromosome (only need austosomes) and apply MAF filter
for chr in {1..22}; do 
	bcftools view -r chr${chr} --min-af 0.05 -Oz -o ${outdir}/gtex_phased_chr${chr}_MAF05.vcf.gz ${gtex_vcf} &
done
wait
# apply MAF filter to 1000 Genomes
for chr in {1..22}; do 
	bcftools view --min-af 0.05 -Oz -o ${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz ${tg}/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz &
done
wait

# change 1000 Genomes naming convention 
for chr in {1..22}; do 
	echo "${chr} chr${chr}" >> chr_name_conv.txt
done
for chr in {1..22}; do 
	bcftools annotate --rename-chrs chr_name_conv.txt -Oz -o ${outdir}/tg_phased_chr${chr}_MAF05_renamed.vcf.gz ${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz &
done

# index all VCFs
rm ${outdir}/tg_phased_chr{1..22}_MAF05.vcf.gz
for file in `ls | grep "gz"`; do 
	tabix -p vcf $file &
done

# find intersection; merge result
for chr in {1..22}; do 
	mkdir -p ${outdir}/chr${chr}
done
wait
for chr in {1..22}; do 
	bcftools isec --collapse none -Oz -p ${outdir}/chr${chr} ${outdir}/gtex_phased_chr${chr}.vcf.gz ${outdir}/tg_phased_chr${chr}_MAF05_renamed.vcf.gz & 
done
wait
for chr in {1..22}; do
	bcftools merge -0 -m none -Oz -o ${outdir}/merged_hg38_chr${chr}.vcf.gz ${outdir}/chr${chr}/0002.vcf.gz ${outdir}/chr${chr}/0003.vcf.gz &
done 
wait


# break point -- merged VCFs from two projects


# there should no longer be any multiallelic variants

################################################################################################################
## Format VCFs for RFMix
################################################################################################################

# filter subjects and re-apply MAF filter 

indir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged
outdir=${indir}/rfmix_maf10_20190814

gtex=${outdir}/gtex-admixed0.9.txt
tg=/users/nicolerg/gtex-admix/metadata/1kg-subject-popcodes.tsv

grep -E "JPT|CHB|CEU|YRI|GWD|MSL|ESN" ${tg} > ${outdir}/tg-refpops.txt

cd ${outdir} 

# remove intermediate files 
rm gtex_phased_chr* tg_phased_chr*
rm -rf chr{1..22}

cat gtex-admixed0.9.txt > rfmix_keep.txt
cut -f1 tg-refpops.txt >> rfmix_keep.txt

# select 1000 Genomes reference panels and admixed GTEx individuals; re-apply MAF > 0.10 filter 
for chr in {1..22}; do 
	bcftools view --force-samples -S rfmix_keep.txt -Oz --min-af 0.1 -o merged_hg38_chr${chr}_filt_MAF10.vcf.gz ../merged_hg38_chr${chr}.vcf.gz &
done 
wait

# remove SNPs with any missing values 
for chr in {1..22}; do 
	bcftools view -g ^miss -Oz -o ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing.vcf.gz ${outdir}/merged_hg38_chr${chr}_filt_MAF10.vcf.gz &
done
wait

# ~5M variants genome-wide

# generate classes, alleles, and snp_position files for RFMix 
for chr in {1..22}; do 
	Rscript ~/gtex-admix/rfmix_to_ancestry_bed/format-rfmix-v2.R ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing.vcf.gz ${tg} ${chr} ${outdir} &
done 
wait

################################################################################################################
## Convert bp to cM
################################################################################################################

outdir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_inputs

for chr in {1..22}; do
	snp_bed=${outdir}/snp_positions-chr${chr}.bed # output from previous step
	mapfile=/mnt/lab_data/montgomery/nicolerg/RESOURCES/genetic-map/hapmap_genetic_maps/hg38/plink.chr${chr}.GRCh38.map # genetic map

	python ~/gtex-admix/rfmix_to_ancestry_bed/interpolate_maps.py ${snp_bed} ${mapfile} ${outdir}/snp-gd-bed-chr${chr}.gz ${outdir}/snp_locations-chr${chr}.gz &
done

wait

################################################################################################################
## Run RFMix
################################################################################################################

outdir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_outputs
indir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_inputs
subj=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/gtex-admixed0.9.txt

mkdir -p ${outdir}/whole-outputs

# this block runs the algorithm separately for each admixed individual
# it iteratively parses the input files to include only one admixed individual and all reference panels

# just do this for the first 10 individuals for now 
for chr in {1..22}; do 
	gunzip ${indir}/snp_locations-chr${chr}.gz
	n_ind=$(cat $subj | wc -l)
	for ind in $(seq 1 $n_ind); do
	#for ind in $(seq 1 10); do
		echo "Working on chromosome ${chr} individual ${ind}..."
		prefix=${outdir}/whole-outputs/gtex-${ind}-chr${chr}
		start=$((${ind}*2-1))
		end=$((${ind}*2))
		tg_start=$(($n_ind*2+1))
		tg_end=$(($(head -1 ${indir}/classes-chr${chr} | wc -w)))
		cut -d' ' -f${start}-${end},${tg_start}-${tg_end} ${indir}/classes-chr${chr} > ${indir}/classes-${chr}-${ind}

		alleles_in=${indir}/alleles-chr${chr}
		cut -c${start}-${end},${tg_start}-${tg_end} ${alleles_in} > ${indir}/alleles-${ind}-chr${chr}
		
		classes=${indir}/classes-${chr}-${ind}
		snp_locations=${indir}/snp_locations-chr${chr}
		alleles=${indir}/alleles-${ind}-chr${chr}

		taskset -c 24-47 python ~/RFMix_v1.5.4/RunRFMix.py \
			PopPhased \
			${alleles} \
			${classes} \
			${snp_locations} \
			'~/RFMix_v1.5.4' \
			--forward-backward \
			--skip-check-input-format \
			-o ${prefix} &

		running=`ps -ef | grep "nicolerg" | grep "RunRFMix.py" | wc -l`
		while [ $running -gt 24 ]; do
			sleep 120
			running=`ps -ef | grep "nicolerg" | grep "RunRFMix.py" | wc -l`
		done 
	done
done















