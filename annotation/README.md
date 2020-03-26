# Annotation

## Calculate variance in gene expression explained by local and global ancestry 
First, use [`interpolate-tss-la.py`](interpolate-tss-la.py) to extract local ancestry calls for 117AX at annotated hg19 TSSs. Then use [`variance_explained_ancestry_gene_expression.R`](variance_explained_ancestry_gene_expression.R) to calculate variance in gene expression explained by each type of ancestry, independent of the other type. This is done using a two-step regression for each gene: first, 117AX gene expression is residualized by one ancestry type (local or global). Then a second linear regression is fit on the residuals, where the predictors are the opposite ancestry type. Variance explained in the expression of a gene is taken to be the R-squared value from the second fit. This process is performed for all expressed genes (all genes in GTEx v8 normalized expression files) for each tissue. 

## Calculate LD between pairs of lead SNPs that are different between LocalAA and GlobalAA  
1. Run [`pairs_for_ld.R`](pairs_for_ld.R) to generate a file `evariant_pairs_test_ld.txt` which includes all of the pairs of SNPs for which to calculate LD  
2. Prepare input files for PLINK: 
  - `cut -f 3,4 evariant_pairs_test_ld.txt | sed "s/  / /g" > pairs_in.txt`  
  - `split -d -n l/15 pairs_in.txt pairs_in_`
3. Run [`plink_ld.sh`](plink_ld.sh) on each `pairs_in_*` chunk to calculate LD between pairs of SNPs using PLINK: 

