#!/bin/bash
# for Stanford SCG (interactive partition)

# generate PEER factors from admixed individuals only

# start with the GTEx normalized expression file 

set -e

tissue=Whole_Blood

EXPRDIR=/path/to/norm/expr # directory with all GTEx v8 normalized expression files 
COVDIR=/path/to/covariate/files
expr=${EXPRDIR}/${tissue}.v8.normalized_expression.bed.gz
outdir=/path/to/outdir

# filter expression down to admixed individuals 
module load r/3.6
Rscript filter_expression_admixed.R ${expr} ../gtex-admix0.9.txt ${outdir}

num_peer=15

mkdir -p ${outdir}/peer

# an old version of R is needed to compile PEER:
module load legacy
module load r/3.3.1
# run_PEER.R is from GTEx eQTL calling pipeline GitHub: 
# https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R
Rscript run_PEER.R \
	${outdir}/${tissue}.v8.normalized_expression.admixed_subset.bed \
	${tissue} \
	${num_peer} \
	-o ${outdir}/peer
# This will generate 3 files:
# ${prefix}_PEER_residuals.txt
# ${prefix}_PEER_alpha.txt
# ${prefix}_PEER_covariates.txt

module purge
module load r/3.6

# make all_covariates files

# add genotype PCs, sex, platform, and pcr to covariate files
# reformat covariate files
output=${outdir}/other_covariates.txt

# GTEx v8 covariates files for each tissue 
cov=${COVDIR}/${tissue}.v8.covariates.txt

head -6 $cov > $output
tail -3 $cov >> $output

Rscript concat_cov.R \
	$output \
	$outdir/peer/${tissue}.PEER_covariates.txt \
	$outdir/${tissue}.all_covariates.txt
