#!/bin/bash

#SBATCH -A account
#SBATCH -c 1
#SBATCH -o slurm/%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mem-per-cpu=10G
#SBATCH -t 10:00:00

set -e 

module load r/3.5
Rscript eqtl_localaa_globalaa.R \
    ${chr} \
    ${tissue} \
    ${exprfile} \
    ${globalcov} \
    ${geno} \
    ${localcov} \
    ${out} 

