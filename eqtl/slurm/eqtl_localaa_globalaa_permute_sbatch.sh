#!/bin/bash

#SBATCH -A smontgom
#SBATCH -p batch
#SBATCH -c 1
#SBATCH -o /oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/slurm/%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=30G
#SBATCH -t 10-00:00:00

set -e 

module load r/3.5
Rscript /oak/stanford/groups/smontgom/nicolerg/src/gtex-admixture-la/eqtl/eqtl_localaa_globalaa_permute.R \
    ${chr} \
    ${tissue} \
    ${exprfile} \
    ${globalcov} \
    ${geno} \
    ${localcov} \
    ${out} \
    ${egenes_master} \
    ${subjid}

