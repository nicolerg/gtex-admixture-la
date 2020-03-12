import gzip

hg19_as_snps = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/ng.3432-S7'
dbsnp_hg38 = '/mnt/lab_data/montgomery/nicolerg/dbsnp/snp150.hg38.txt.gz'
outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/ng.3432-S7_hg38.txt'

allele_sp_hg19 = {}

with open(hg19_as_snps) as hg19:
	next(hg19)
	for line in hg19:
		l = line.strip().split()
		rsid = l[3].split(',')[0]
		pred.fit.pctSig = l[4]
		genotypes = l[8]
