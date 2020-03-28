#!/bin/bash
# variant_id    panel_variant_id        chromosome      position        effect_allele   non_effect_allele    current_build    frequency       sample_size     zscore  pvalue  effect_size     standard_error  imputation_status     n_cases

set -e

# cd /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_gwas
# ls imputed_gwas_hg38_1.1 | grep -E "^imputed" | parallel --jobs 20 "bash /users/nicolerg/gtex-admixture-la/colocalization/format_gtex_gwas.sh {1}" &
file=$1

indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_gwas
cd ${indir}
mkdir -p coloc_input

# variant_id      gene    chr     snp_pos ref     effect_allele   build   pval    beta    se      effect_af

# imputed_UKB_6152_8_diagnosed_by_doctor_Asthma.txt.gz
new_file=$(echo $file | sed "s/\.txt\.gz/\.formatted\.txt/")
echo "chr	snp_pos	alt	ref	beta	se	pvalue	effect_af" > coloc_input/${new_file}
zcat imputed_gwas_hg38_1.1/${file} | sed -e '1d' | awk '{print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $12 "\t" $13 "\t" $11 "\t" $8}' | sort -k 1,1 -k2,2n >> coloc_input/${new_file}
bgzip -f coloc_input/$new_file
tabix -f -S 1 -s 1 -b 2 -e 2 coloc_input/${new_file}.gz
