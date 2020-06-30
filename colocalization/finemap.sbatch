#!/bin/bash

#SBATCH -p interactive
#SBATCH -c 10
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=100G
#SBATCH -t 10-00:00:00
#SBATCH -J adipose_global

set -e

## to submit all jobs:
# cd config/finemap/
# for config_prefix in $(ls *finemap.json | sed "s/_finemap.*//"); do
# 	sed "s/JOBNAME/${config_prefix}/" ../../finemap.sbatch > ${config_prefix}_finemap.sbatch
# 	config=${config_prefix}_finemap.json
# 	sbatch --export=config=${config} ${config_prefix}_finemap.sbatch
# done

module load r/3.6 # R v3.6
module load miniconda/2 # python2
module load plink/1.90b6.13 # PLINK v1.90b
module load tabix # tabix

# global variable: config=config/finemap/adipose_global_finemap.json
python dispatch.py ${config} 10
