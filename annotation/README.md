# Annotation

## Calculate variance in gene expression explained by local and global ancestry 
1. Use [`interpolate-tss-la.py`](interpolate-tss-la.py) to extract local ancestry calls for 117AX at annotated hg19 TSSs 
2. Use [`variance_explained_ancestry_gene_expression.R`](variance_explained_ancestry_gene_expression.R) to calculate variance in gene expression explained by each type of ancestry, independent of the other type. This is done using a two-step regression for each expressed gene in each tissue (all genes in GTEx v8 normalized expression files):
  1. First, 117AX gene expression for a single gene is residualized by one ancestry type (local or global)
  2. A second linear regression is fit on the residuals, where the predictors are the opposite ancestry type. Variance explained in the expression of a gene is taken to be the R-squared value from the second fit.  

## Calculate LD between pairs of lead SNPs that are different between LocalAA and GlobalAA  
1. Run [`pairs_for_ld.R`](pairs_for_ld.R) to generate a file `evariant_pairs_test_ld.txt` which includes all of the pairs of SNPs for which to calculate LD  
2. Prepare input files for PLINK: 
    ```bash 
    cut -f 3,4 evariant_pairs_test_ld.txt | sed "s/ / /g" > pairs_in.txt
    split -d -n l/15 pairs_in.txt pairs_in_ 
    ```
3. Run [`plink_ld.sh`](plink_ld.sh) on each `pairs_in_*` chunk to calculate LD between pairs of SNPs using PLINK: `for file in pairs_in_*; do bash plink_ld.sh ${file} & done`
4. Concatenate results  
5. Run [`define_diff_lead_by_ld.R`](define_diff_lead_by_ld.R) to restrict the definition of "different lead SNPs" to those with LD < 1  

## Allele frequencies from 1000 Genomes 
See [`get_1kg_af.sh`](get_1kg_af.sh), which calls [`match_af_gtex_id.py`](match_af_gtex_id.py).  

## Map GTEx variant IDs to RSIDs
There are two ways to do this:  

### Use dbsnp:
1. Download the dbsnp database for [hg19](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp150.txt.gz) and/or [hg38](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/snp150.txt.gz)  
2. Use [`gtex_v8_rsid_map.py`](gtex_v8_rsid_map.py) to generate a map between hg19 RSIDs, hg38 RSIDs, and hg38 GTEx variant IDs 

### Use 1000 Genomes VCFs (hg38): 
See [`rsid_from_1kg_vcf.py`](rsid_from_1kg_vcf.py) for one way to map between RSIDs used as 1000 Genomes variant IDs and GTEx v8 variant IDs. Alternatively, use `bcftools` to merge GTEx and 1000 Genomes VCFs [here](https://github.com/nicolerg/gtex-admixture-la/tree/master/rfmix#merge-gtex-and-1000-genomes-vcfs) and pull the RSID-to-GTEx variant ID from the `ID` column, e.g. `chr13_18445955_C_T_b38;rs140871821`.  

## Variance in SNP genotype explained by local ancestry 
Use [`r2_all_snps.R`](r2_all_snps.R) to get variance in genotype explained by local ancestry for all lead SNPs in the LocalAA/GlobalAA eQTL analysis (117AX only).  

Then use [`gtex_variant_la_r2_by_tissue.R`](gtex_variant_la_r2_by_tissue.R) to calculate variance in genotype explained by local ancestry for all GTEx v8 eSNPs (all subjects).  
