#!/bin/bash
# LiftoverVcf for GTEx: hg38 --> hg19

refhg19=/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa

input=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
outdir=/mnt/lab_data/montgomery/nicolerg/filtered
output=${outdir}/gtex.phased.MAF01.hg19.vcf
reject=${outdir}/rejected_variants.vcf
chain38to19=/mnt/lab_data/montgomery/shared/liftOver/chains/hg38ToHg19.over.chain.gz
zcat ${chain38to19} > ${outdir}/hg19ToHg38.over.chain

taskset -c 24-35 java -jar ~/picard.jar LiftoverVcf \
     I=${input} \
     O=${output} \
     CHAIN=${outdir}/hg19ToHg38.over.chain \
     REJECT= ${reject}\
     R=${refhg19} # target reference build

