#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm/slurm-%j.out

set -e 

module load r/3.6 # R v3.6
module load miniconda/2 # python2
module load plink/1.90b6.13 # PLINK v1.90b
module load tabix # tabix

config=/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/overlap/finemap_1000G_sig.json
python brain_gwas/scripts/dispatch.py ${config} 12
