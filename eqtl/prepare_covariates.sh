#!/bin/bash

# generate PEER factors from admixed individuals only

# start with the GTEx normalized expression file 

tissue=Skin_Not_Sun_Exposed_Suprapubic

# EXPRDIR: directory with all GTEx v8 normalized expression files 
expr=${EXPRDIR}/${tissue}.v8.normalized_expression.bed.gz
outdir=/path/to/outdir

# filter expression down to admixed individuals 
Rscript filter_expression_admixed.R ${expr} ${outdir}

num_peer=15

mkdir -p ${outdir}/peer

# run_PEER.R is from GTEx eQTL calling pipeline GitHub: 
# https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R
taskset -c 12-47 Rscript run_PEER.R \
	${outdir}/${tissue}.v8.normalized_expression.admixed_subset.bed \
	${tissue} \
	${num_peer} \
	-o ${outdir}/peer
# This will generate 3 files:
# ${prefix}_PEER_residuals.txt
# ${prefix}_PEER_alpha.txt
# ${prefix}_PEER_covariates.txt

# make all_covariates files

# add genotype PCs, sex, platform, and pcr to covariate files
# reformat covariate files
output=${outdir}/other_covariates.txt

# GTEx v8 covariates files for each tissue 
cov=GTEx_Analysis_v8_eQTL_covariates/${tissue}.v8.covariates.txt

head -6 $cov > $output
tail -3 $cov >> $output

Rscript concat_cov.R \
	$output \
	$outdir/peer/${tissue}.PEER_covariates.txt \
	$outdir/${tissue}.all_covariates.txt
