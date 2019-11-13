## eQTL calling 

### Prepare VCFs



filter_vcf.sh applies the MAC10 filter 

prepare_covariates.sh
calls:
  - filter_expression_admixed.R
  - concat_cov.R

batch_eqtl_localaa_globalaa.sh
calls:
  - eqtl_localaa_globalaa.R
