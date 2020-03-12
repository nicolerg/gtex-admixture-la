# Print a table of the RACE codes for the individuals in the GTEx v8 vcf (GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz)

## Import phenotype tables
# GTEx 
# RACE column:
#	1 = Asian
# 	2 = Black or African American
# 	3 = White
# 	4 = American Indian or Alaska Native
# 	99 = Unknown
# 	98 = Not Reported
phenofile <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
gpheno <- read.delim(file=phenofile, sep='\t', header=TRUE, stringsAsFactors=FALSE)
gtexfile <- "/mnt/lab_data/montgomery/nicolerg/filtered/gtex.chr21.vcf.gz"
gtex <- read.delim(gtexfile, sep='\t', header=FALSE, nrows=1, stringsAsFactors=FALSE)
names <- gtex[1,]
m <- match(names, gpheno$SUBJID)
table(gpheno$RACE[m])