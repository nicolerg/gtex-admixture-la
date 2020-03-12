#!/bin/bash 

############################################################################################################################
## get pairwise Fst values between EUR subpops and AFR subpops
############################################################################################################################ 

### adjusted to test just RSIDs that were missing in the first round

outdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/round2
srcdir=/users/nicolerg/gtex-admix/scripts/LAVA

snp_list=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/need_to_test_rsid.txt

for i in $(seq 1 $(cat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/pop1.txt | wc -l)); do 
	pop1=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/$(awk 'NR == a' a=$i /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/pop1.txt)
	pop2=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/$(awk 'NR == a' a=$i /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/pop2.txt)

	pop1_short=$(basename $pop1 | sed "s/_subjs.*//")
	pop2_short=$(basename $pop2 | sed "s/_subjs.*//")

	for chr in {1..22}; do 
		taskset -c 0-21 vcftools --gzvcf /mnt/lab_data/montgomery/shared/1KG/hg38/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz --snps ${snp_list} --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out ${outdir}/${pop1_short}-${pop2_short}-chr${chr} &
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
	cat ${outdir}/${i}-chr{1..22}.weir.fst >> ${outdir}/${i}.weir.fst
done
rm ${outdir}/*-chr{1..22}.weir.fst
mkdir -p ${outdir}/log 
mv ${outdir}/*.log ${outdir}/log
Rscript ${srcdir}/annotation/compile_fst.R 

# merge new and old results --> /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master-merged_eur_afr_all_snps.RData

# combine Fst and egenes
Rscript ${srcdir}/individual_final_figs.R --egene_fst --cutoff=1e-06
Rscript ${srcdir}/individual_final_figs.R --egene_fst --cutoff=1e-05
Rscript ${srcdir}/individual_final_figs.R --egene_fst --cutoff=1e-04
