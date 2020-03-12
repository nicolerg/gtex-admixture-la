#!/bin/bash

# pipeline for analyzing eQTL calling results, starting from merged allpairs file 
# 22 May 2019
# Nicole Gay

set -e 

srcdir=/users/nicolerg/gtex-admix/scripts/LAVA
outdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed

############################################################################################################################
## get all lead SNPs; filter allpairs file down to tests with pval < 0.05 in either method for ease of downstream analysis
############################################################################################################################

for tissue in Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal Nerve_Tibial Artery_Tibial Lung Adipose_Subcutaneous; do

	allpairs=${outdir}/${tissue}/${tissue}-LAVA-global-allpairs-merged.tsv.gz

	if [ -f "${outdir}/${tissue}/gtex.admix.global.egenes.tied.txt.gz" ]; then 
		:
	else
		# get lead SNPs for each method
		python ${srcdir}/extract-egene-lava.py ${allpairs} ${tissue} ${outdir}/${tissue} &
		# ${outdir}/${tissue}/gtex.admix.global.egenes.tied.txt.gz
		# ${outdir}/${tissue}/gtex.admix.lava.egenes.tied.txt.gz
	fi

	if [ -f "${outdir}/${tissue}/${tissue}-LAVA-global-allpairs-filt-merged.tsv.gz" ]; then 
		:
	else
		# filter allpairs file by p-value 
		python ${srcdir}/filter_allpairs.py ${allpairs} &
		# ${outdir}/${tissue}/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-allpairs-filt-merged.tsv.gz
	fi

done

wait 

############################################################################################################################
## now, make egenes_master and overlap files with all tissues in merged files (this is hard-coded for the tissues I want)
## these are sort of the "master" files used for all of the annotation scripts
############################################################################################################################

if [ ! -f "${outdir}/merged/same_egene_diff_lead_1e-04.RData" ]; then 
	Rscript ${srcdir}/annotation/get_uniq_egenes.R ${outdir} &
fi

############################################################################################################################
## get a quick summary of eGenes results (counts and overlap) at a given p-value cutoff 
############################################################################################################################

# eQTL calling summary 
Rscript ${srcdir}/annotation/eqtl_calling_summary.R 1e-6
Rscript ${srcdir}/annotation/eqtl_calling_summary.R 1e-5

############################################################################################################################
## make all of the plots in figure 2:
# q-q plot
# egene discovery per tissue
# egene replication distribution
# overlap of eGenes
############################################################################################################################

Rscript ${srcdir}/individual_final_figs.R --qq & 
Rscript ${srcdir}/individual_final_figs.R --tss_distance --egene_discovery --pval_distn --venn --cutoff=1e-5 &
Rscript ${srcdir}/individual_final_figs.R --tss_distance --egene_discovery --pval_distn --venn --cutoff=1e-6 &

############################################################################################################################
## ASB (plot and empiral p-values)
############################################################################################################################

if [ ! -f "$outdir/annotation/allele-specific-binding/gtex_mac10_snp_to_asb_map.tsv" ]; then 
	python ${srcdir}/annotation/allele_specific_binding.py $outdir/annotation 
fi
Rscript ${srcdir}/individual_final_figs.R --asb --cutoff=1e-05
Rscript ${srcdir}/individual_final_figs.R --asb --cutoff=1e-06
Rscript ${srcdir}/individual_final_figs.R --asb --cutoff=1e-04
# permutation test for empirical p-value 
for tissue in Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal Nerve_Tibial Artery_Tibial Lung Adipose_Subcutaneous; do
	python ${srcdir}/annotation/asb_perm.py ${outdir}/annotation 1e-06 $tissue &
	python ${srcdir}/annotation/asb_perm.py ${outdir}/annotation 1e-05 $tissue & 
done
wait
# manually edit p-value labels in individual_final_figs.R 
############################################################################################################################
## TF binding 
############################################################################################################################

for tissue in Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal Nerve_Tibial Artery_Tibial Lung Adipose_Subcutaneous; do 
	python ${srcdir}/annotation/encode_tf_binding.py $tissue &
done
wait
Rscript ~/gtex-admix/scripts/LAVA/individual_final_figs.R --tf_bind --cutoff=1e-06
Rscript ~/gtex-admix/scripts/LAVA/individual_final_figs.R --tf_bind --cutoff=1e-05

############################################################################################################################
## CADD
############################################################################################################################ 

python ${srcdir}/annotation/parse_cadd.py 
Rscript ~/gtex-admix/scripts/LAVA/individual_final_figs.R --cadd_plot --cutoff=1e-06 # only Tibial Nerve is significantly different
Rscript ~/gtex-admix/scripts/LAVA/individual_final_figs.R --cadd_plot --cutoff=1e-05 # no tissues are significantly different 

############################################################################################################################
## correlation between SNP genotypes and local ancestry; stratify by differences in population-level AFs 
############################################################################################################################ 

# get R-squared values between all lead SNPs (both my version and GTEx v8) and local ancestry 
for chr in {1..22}; do 
	Rscript ${srcdir}/annotation/r2_all_snps.R $chr 1 &
done 
wait
# get allele frequencies
if [ ! -f  "$outdir/annotation/af/gtex.admixed.MAC10.snps.af.txt.gz" ]; then 
	mkdir -p ${outdir}/annotation/af
	bash ${srcdir}/annotation/get_1kg_af.sh ${outdir}/annotation/af
fi
# merge R2
if [ ! -f "$outdir/annotation/r2/all_lead_snp_r2.txt" ]; then 
	head -1 ${outdir}/annotation/r2/1.r2.txt > ${outdir}/annotation/r2/all_lead_snp_r2.txt
	sed -e '1d' ${outdir}/annotation/r2/{1..22}.r2.txt >> ${outdir}/annotation/r2/all_lead_snp_r2.txt
	gzip ${outdir}/annotation/r2/all_lead_snp_r2.txt
fi
# merge all r2-af files; plot 
Rscript ${srcdir}/annotation/r2_all_snps.R 1
Rscript ${srcdir}/individual_final_figs.R --local_var &

############################################################################################################################
## get correlation between GTEx leads SNPs and local ancestry, ALL INDIVIDUALS!! 
############################################################################################################################ 

# NOT DONE 
# write a new script for this. borrow from ${srcdir}/annotation/r2_all_snps.R 

############################################################################################################################
## get pairwise Fst values between EUR subpops and AFR subpops
############################################################################################################################ 
for subpop in $(grep "	AFR" /users/nicolerg/gtex-admix/metadata/1kg-subject-popcodes.tsv | cut -f2 | sed -e '1d' | sort | uniq); do 
	grep "	$subpop	" /users/nicolerg/gtex-admix/metadata/1kg-subject-popcodes.tsv | cut -f1 > ${outdir}/annotation/fst/AFR_${subpop}_subjs.txt
done

for subpop in $(grep "	EUR" /users/nicolerg/gtex-admix/metadata/1kg-subject-popcodes.tsv | cut -f2 | sed -e '1d' | sort | uniq); do 
	grep "	$subpop	" /users/nicolerg/gtex-admix/metadata/1kg-subject-popcodes.tsv | cut -f1 > ${outdir}/annotation/fst/EUR_${subpop}_subjs.txt
done

# # manual list of pairs 
# pop1.txt:	        pop2.txt:
# AFR_ESN_subjs.txt AFR_ACB_subjs.txt
# AFR_YRI_subjs.txt AFR_MSL_subjs.txt
# AFR_YRI_subjs.txt AFR_LWK_subjs.txt
# AFR_YRI_subjs.txt AFR_GWD_subjs.txt
# AFR_YRI_subjs.txt AFR_ASW_subjs.txt
# AFR_YRI_subjs.txt AFR_ESN_subjs.txt
# AFR_YRI_subjs.txt AFR_ACB_subjs.txt
# AFR_MSL_subjs.txt AFR_LWK_subjs.txt
# AFR_MSL_subjs.txt AFR_GWD_subjs.txt
# AFR_MSL_subjs.txt AFR_ASW_subjs.txt
# AFR_MSL_subjs.txt AFR_ESN_subjs.txt
# AFR_MSL_subjs.txt AFR_ACB_subjs.txt
# AFR_LWK_subjs.txt AFR_GWD_subjs.txt
# AFR_LWK_subjs.txt AFR_ASW_subjs.txt 
# AFR_LWK_subjs.txt AFR_ESN_subjs.txt
# AFR_LWK_subjs.txt AFR_ACB_subjs.txt
# AFR_GWD_subjs.txt AFR_ASW_subjs.txt
# AFR_GWD_subjs.txt AFR_ESN_subjs.txt
# AFR_GWD_subjs.txt AFR_ACB_subjs.txt
# AFR_ASW_subjs.txt AFR_ESN_subjs.txt
# AFR_ASW_subjs.txt AFR_ACB_subjs.txt
# EUR_FIN_subjs.txt EUR_CEU_subjs.txt
# EUR_TSI_subjs.txt EUR_IBS_subjs.txt
# EUR_TSI_subjs.txt EUR_GBR_subjs.txt
# EUR_TSI_subjs.txt EUR_FIN_subjs.txt
# EUR_TSI_subjs.txt EUR_CEU_subjs.txt
# EUR_IBS_subjs.txt EUR_GBR_subjs.txt
# EUR_IBS_subjs.txt EUR_FIN_subjs.txt
# EUR_IBS_subjs.txt EUR_CEU_subjs.txt
# EUR_GBR_subjs.txt EUR_FIN_subjs.txt
# EUR_GBR_subjs.txt EUR_CEU_subjs.txt

snp_list=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.rsid.txt
zcat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.rsid.map.txt.gz | sed -e '1d' | cut -f1 > ${snp_list}
## NOTE: I used an older snp_list to run this pipeline and ended up needing to exclude all results where there were duplicate positions
##		EUR: 6148995 1000 Genomes SNPs (excluding NA and duplicate positions)
##		AFR: 7350438 1000 Genomes SNPs (excluding NA and duplicate positions)
## the current RSID map includes the correct RSID by position according to the alternate allele
## there are RSIDs for 7550023 out of 8088666 SNPs tested for admixed individuals 
for i in $(seq 1 $(cat ${outdir}/annotation/fst/pop1.txt | wc -l)); do 
	pop1=${outdir}/annotation/fst/$(awk 'NR == a' a=$i ${outdir}/annotation/fst/pop1.txt)
	pop2=${outdir}/annotation/fst/$(awk 'NR == a' a=$i ${outdir}/annotation/fst/pop2.txt)

	pop1_short=$(basename $pop1 | sed "s/_subjs.*//")
	pop2_short=$(basename $pop2 | sed "s/_subjs.*//")

	for chr in {1..22}; do 
		taskset -c 0-23 vcftools --gzvcf /mnt/lab_data/montgomery/shared/1KG/hg38/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz --snps ${snp_list} --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out ${outdir}/annotation/fst/${pop1_short}-${pop2_short}-chr${chr} &
	done 
	wait 
done
wait

comparisons=( AFR_ASW-AFR_ACB \
	AFR_ASW-AFR_ESN \
	AFR_ESN-AFR_ACB \
	AFR_GWD-AFR_ACB \
	AFR_GWD-AFR_ASW \
	AFR_GWD-AFR_ESN \
	AFR_LWK-AFR_ACB \
	AFR_LWK-AFR_ASW \
	AFR_LWK-AFR_ESN \
	AFR_LWK-AFR_GWD \
	AFR_MSL-AFR_ACB \
	AFR_MSL-AFR_ASW \
	AFR_MSL-AFR_ESN \
	AFR_MSL-AFR_GWD \
	AFR_MSL-AFR_LWK \
	AFR_YRI-AFR_ACB \
	AFR_YRI-AFR_ASW \
	AFR_YRI-AFR_ESN \
	AFR_YRI-AFR_GWD \
	AFR_YRI-AFR_LWK \
	AFR_YRI-AFR_MSL \
	EUR_FIN-EUR_CEU \
	EUR_GBR-EUR_CEU \
	EUR_GBR-EUR_FIN \
	EUR_IBS-EUR_CEU \
	EUR_IBS-EUR_FIN \
	EUR_IBS-EUR_GBR \
	EUR_TSI-EUR_CEU \
	EUR_TSI-EUR_FIN \
	EUR_TSI-EUR_GBR \
	EUR_TSI-EUR_IBS )

for i in "${comparisons[@]}"; do
	cat ${outdir}/annotation/fst/${i}-chr{1..22}.weir.fst >> ${outdir}/annotation/fst/${i}.weir.fst
done
rm ${outdir}/annotation/fst/*-chr{1..22}.weir.fst
mkdir -p ${outdir}/annotation/fst/log 
mv ${outdir}/annotation/fst/*.log ${outdir}/annotation/fst/log
Rscript ${srcdir}/annotation/compile_fst.R 

# combine Fst and egenes
Rscript ${srcdir}/individual_final_figs.R --egene_fst --cutoff=1e-06
Rscript ${srcdir}/individual_final_figs.R --egene_fst --cutoff=1e-05
Rscript ${srcdir}/individual_final_figs.R --egene_fst --cutoff=1e-04

############################################################################################################################
## colocalization for all genes where lead SNP is different (1e-4)
############################################################################################################################ 

bash ${srcdir}/coloc/coloc_pipeline.sh 
Rscript ${srcdir}/individual_final_figs.R --coloc

# combine coloc and fst 
Rscript ${srcdir}/individual_final_figs.R --coloc_fst 



