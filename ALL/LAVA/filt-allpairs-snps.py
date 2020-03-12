import sys
import gzip 

# filter allpairs files by SNPs in filtered VCFs

vcfdir = '/mnt/lab_data/montgomery/nicolerg/genotypes-hg19/gtex-admixed-hg19'
pairsdir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/by-chr'
outdir = pairsdir + '/maf-filt'

chr = sys.argv[1]
invcf = vcfdir + '/gtex.admixed.maf.chr' + chr + '.recode.vcf.gz'
inpairs = pairsdir + '/lung-LAVA-global-merged-allpairs-chr' + chr + '.tsv.gz'
outfile = outdir + '/lung-LAVA-merged-maffilt-allpairs-chr' + chr + '.tsv.gz'

snpset = set()

# get the set of filtered SNPs
with gzip.open(invcf, 'rb') as vcf:
	next(vcf)
	for line in vcf:
		l = line.strip().split()
		snp = l[2]
		snpset.add(snp)

header = '\t'.join(['gene_id', \
	'variant_id', \
	'tss_distance', \
	'ma_samples', \
	'ma_count', \
	'global_pval_nominal', \
	'global_slope', \
	'global_slope_se', \
	'lava_pval_nominal', \
	'lava_slope', \
	'lava_slope_se', \
	'afr_pval_nominal', \
	'afr_slope', \
	'afr_slope_se']) + '\n'

# now filter allpairs file
with gzip.open(inpairs, 'rb') as inp, gzip.open(outfile, 'wb') as out:
	next(inp)
	out.write(header)
	for line in inp:
		line = line.strip()
		snp = line.split()[1]
		if snp in snpset:
			out.write(line+'\n')
