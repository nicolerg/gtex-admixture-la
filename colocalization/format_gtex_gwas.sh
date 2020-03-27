#!/bin/bash
# variant_id	panel_variant_id	chromosome	position	effect_allele	non_effect_allele	current_build	frequency	sample_size	zscore	pvalue	effect_size	standard_error	imputation_status	n_cases

set -e

# variant_id      gene    chr     snp_pos ref     effect_allele   build   pval    beta    se      effect_af
for file in `ls | grep -E "^imputed"`; do 

	# imputed_UKB_6152_8_diagnosed_by_doctor_Asthma.txt.gz
	new_file=$(echo $file | sed -e "s/imputed/coloc_imputed/" -e "s/\.gz//") 
	echo "chr	snp_pos	alt	ref	beta	se	pvalue	effect_af" > coloc_input/${new_file}
	zcat $file | sed -e '1d' | awk '{print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $12 "\t" $13 "\t" $11 "\t" $8}' | sort -k 1,1 -k2,2n >> coloc_input/${new_file}
	bgzip -f coloc_input/$new_file
	tabix -f -S 1 -s 1 -b 2 -e 2 coloc_input/${new_file}.gz

done
 
