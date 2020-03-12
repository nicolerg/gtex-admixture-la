# eQTL calling with local ancestry

## Calculate new peer factors for admixed individuals only and recompile global covariates file
* `/users/nicolerg/gtex-admix/scripts/LAVA/eqtl-calling/new_peer_admixed.sh`
  Calls:
  * `/users/nicolerg/gtex-admix/scripts/gtex-qtl-pipeline/src/run_PEER.R`
  * `/users/nicolerg/gtex-admix/scripts/LAVA/eqtl-calling/concat_cov.R`
  * `/users/nicolerg/gtex-admix/scripts/LAVA/eqtl-calling/filter_expression_admixed.R`

## Run eQTL calling 
Root directory: `~/gtex-admix/scripts/LAVA/eqtl-calling`
| File 					| Description					|
|-----------------------|-------------------------------|
| `gtex-lava-global.R`  | Most updated Rscript for standard pipeline |
| `batch-master-lava-global.sh` | Most updated shell script to run eQTL calling pipeline |
| `batch-master-local-eqtl-allsubj.sh` | SCG `sbatch`; run eQTL all calling with local ancestry and all 838 subjects |
| `batch-master-local-eqtl.sh` | SCG `sbatch`; run eQTL calling with both local ancestry and global genotype PCs on just the 117 admixed invididuals |
| `global-local-eqtl-custom-peer.sh` | Run `gtex-lava-global.sh` with new PEER factors calculated only with 117 admixed subjects |
| `gtex-lava-global-allsubj.R` | Called by `batch-master-local-eqtl-allsubj.sh` |
| `gtex-lava-global-both.R` | Called by `gtex-lava-global-both-scg.sh` |
| `gtex-lava-global-both-scg.sh` | Called by `batch-master-local-eqtl.sh` |

## ID different sets of unique and overlapping eGenes between Local and Global methods
1. Filter down allpairs files: `/users/nicolerg/gtex-admix/scripts/LAVA/filter_allpairs.py`
2. Get lead SNPs per eGene from the original allpairs file: `/users/nicolerg/gtex-scripts/LAVA/extract-egene-lava.py`
3. Run `/users/nicolerg/gtex-scripts/LAVA/annotation/get_uniq_egenes.R` with the correct file paths 

## Chromosome 19 comparisons
1. Make merged files with p-value, slope, and count from each iteration of eQTL calling: `/users/nicolerg/gtex-admix/scripts/LAVA/chr19_combined_egene.py`
2. Make plots of the overlap: `/users/nicolerg/gtex-admix/scripts/LAVA/chr19_comparisons.R`

## Whole genome
1. Get lead SNPs per eGene: `/users/nicolerg/gtex-scripts/LAVA/extract-egene-lava.py`
2. Make egenes_master.RData 
  i. Filter down allpairs to pval < 0.05: `/users/nicolerg/gtex-admix/scripts/LAVA/filter_allpairs.py`
3. Compare power between the two methods 
4. Look at overlap in eQTLs
5. TF binding
6. ASB 
7. TSS distance
8. Colocalization 
9. Correlation between genotype of lead SNPs and local ancestry
10. HOMER (use DeepBind or DeepSEA instead?)
11. 