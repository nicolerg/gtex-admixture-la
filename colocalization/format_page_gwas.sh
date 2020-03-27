#!/bin/bash

set -e

# Chr 
# Position_hg19 
# SNP 
# Other-allele 
# Effect-allele 
# Effect-allele-frequency 
# Sample-size 
# Effect-allele-frequency-cases 
# Sample-size-cases 
# Beta 
# SE 
# P-val 
# INFO-score 
# rsid

## to run:
# cd /mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/coloc/PAGE
# ls WojcikG_PMID_* | parallel --jobs 12 "bash /users/nicolerg/gtex-admixture-la/colocalization/format_page_gwas.sh {1}" &

file=$1

indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/REVISIONS/coloc/PAGE
cd ${indir}

mkdir -p bed liftover coloc_input 

# needs to be converted from hg19 to hg38 using LiftOver
# convert to BED 
new_file=$(echo $file | sed -e "s/\.gz/\.bed/") 
#echo "chr	snp_pos_1	snp_pos	alt;ref;beta;se;pvalue;effect_af" > bed/${new_file}
zcat $file | sed -e '1d' | awk '{printf "chr%s\t%.0f\t%.0f\t%s;%s;%.7f;%.8f;%.8f;%.7f\n",$1,$2-1,$2,$5,$4,$10,$11,$12,$6}' > bed/${new_file}

# convert to hg38 
hg38=$(echo ${new_file} | sed "s/bed/hg38\.bed/")
liftOver bed/${new_file} liftover/hg19ToHg38.over.chain.gz liftover/${hg38} liftover/${hg38}_unMapped

## reformat for coloc input 
formatted=$(echo ${hg38} | sed "s/hg38/hg38\.formatted/")
cut -f 1,3,4 liftover/${hg38} >> coloc_input/${formatted}
sed -i $"s/;/\t/g" coloc_input/${formatted}
echo -e "chr\tsnp_pos\talt\tref\tbeta\tse\tpvalue\teffect_af_gwas" > coloc_input/${formatted}.sorted
# sort
sort -k 1,1 -k2,2n coloc_input/${formatted} >> coloc_input/${formatted}.sorted
final=$(echo ${formatted}.sorted | sed "s/\.bed.*/\.sorted\.tsv/")
mv coloc_input/${formatted}.sorted coloc_input/${final}
# bgzip and index
bgzip -f coloc_input/${final}
tabix -f -S 1 -s 1 -b 2 -e 2 coloc_input/${final}.gz
	