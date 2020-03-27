#!/bin/bash

# COLOC eQTL calling pipeline for Global and Local results
# all tissues; same eGene, different lead SNP at 1e-04
# 114 GWAS from coloc preprint + 30 PAGE traits 

srcdir=~/gtex-admixture-la/colocalization

indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl
base=${indir}/REVISIONS/coloc

# list of same egenes, diff lead SNP output by get_uniq_egenes.R 
gene_list=${indir}/REVISIONS/merged/same_egene_diff_lead_1e-04.txt

af=${indir}/admixed/geno/annotation/snp_to_effect_af.tsv.gz
if [ ! -f "${af}" ]; then 
	# get AF of every SNP tested in admixed individuals (AF only within admixed population)
	for chr in {1..22}; do 
		vcf=${indir}/admixed/geno/gtex.admixed.MAC10.phased.chr${chr}.recode.vcf.gz
		vcftools --gzvcf $vcf --out geno/af/${chr}.admixed.af --freq &
	done
	wait 
	# merge with variant_id and write map: geno/snp_to_effect_af.tsv.gz
fi

# parse allpairs with annotation/parse_allpairs_for_coloc.py 
mkdir -p ${base}/allpairs_difflead
for tissue in Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal Nerve_Tibial Artery_Tibial Lung Adipose_Subcutaneous Whole_Blood; do 
	python ${srcdir}/parse_allpairs_for_coloc.py ${tissue} ${gene_list} ${indir}/admixed ${af} ${base}/allpairs_difflead &
done
wait

# sort; bgzip
for file in `ls ${base}/allpairs_difflead | grep "allpairs.tsv"`; do # Skin_Not_Sun_Exposed_Suprapubic.filtered.local.allpairs.tsv.gz
	new_file=$(echo $file | sed "s/allpairs.*/allpairs\.sorted\.tsv/")
	echo ${new_file}
	zcat ${base}/allpairs_difflead/$file | head -1 > ${base}/allpairs_difflead/${new_file}
	zcat ${base}/allpairs_difflead/$file | sed -e '1d' | sort -k 1,1 -k2,2n >> ${base}/allpairs_difflead/${new_file}
	bgzip -f ${base}/allpairs_difflead/${new_file} & # Skin_Not_Sun_Exposed_Suprapubic.filtered.local.allpairs.sorted.tsv.gz
done
wait
for file in `ls ${base}/allpairs_difflead | grep -E "allpairs.sorted.tsv.gz$"`; do
	tabix -f -S 1 -s 1 -b 2 -e 2 ${base}/allpairs_difflead/${file} & 
done
wait

# make config file 
config=${base}/config/117admixed_newpeer.config 

################################################################################################################
# # make coloc/config/n_gwas and coloc/config/n_eqtl sample size files manually

# # format for n_gwas.txt:
# coloc_imputed_ADIPOGen_Adiponectin.txt.gz	29304
# coloc_imputed_Astle_et_al_2016_Eosinophil_counts.txt.gz	173480 

# # format for n_eqtl.txt:
# Adipose_Subcutaneous.filtered.global.allpairs.sorted.tsv.gz	84
# Adipose_Subcutaneous.filtered.local.allpairs.sorted.tsv.gz	84
# Artery_Tibial.filtered.global.allpairs.sorted.tsv.gz	89
# Artery_Tibial.filtered.local.allpairs.sorted.tsv.gz	89
################################################################################################################

# run coloc 
cd ${srcdir}
taskset -c 30-47 python ./dispatch.py ${config} 100 # dispatch.py is one of Mike Gloudeman's scripts 

colocdir=${base}/output/gtex-admixed-difflead/2019-05-24_14-28-31.840186_117admixed
head -1 ${colocdir}/coloc_imputed_UKB_6152_9_diagnosed_by_doctor_Hayfever_allergic_rhinitis_or_eczema_txt_gz_coloc_status.txt > ${colocdir}/merged_coloc.txt
for file in `ls ${colocdir} | grep 'coloc_status.txt'`; do sed -e '1d' ${colocdir}/$file >> ${colocdir}/merged_coloc.txt ; done

