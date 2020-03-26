# Annotation

## Calculate variance in gene expression explained by local and global ancestry 
First, use [`interpolate-tss-la.py`](interpolate-tss-la.py) to extract local ancestry calls for 117AX at annotated hg19 TSSs. Then use [`variance_explained_ancestry_gene_expression.R`](variance_explained_ancestry_gene_expression.R) to calculate variance in gene expression explained by each type of ancestry, independent of the other type. This is done using a two-step regression for each gene: first, 117AX gene expression is residualized by one ancestry type (local or global). Then a second linear regression is fit on the residuals, where the predictors are the opposite ancestry type. Variance explained in the expression of a gene is taken to be the R-squared value from the second fit. This process is performed for all expressed genes (all genes in GTEx v8 normalized expression files) for each tissue. 

