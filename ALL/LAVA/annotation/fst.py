# parse Fst values 

import sys 
import gzip 
import pandas

fst = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/Fst_CEU_vs_YRI.whole_genome.pvalues.gz'
rsid_map = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/rsid_to_gtex_snp_map.tsv'

rsid_snp = {}
gtex_fst = {}

with open(rsid_map, 'rb') as snpmap:
	for line in snpmap:
		l = line.strip().split('\t')
		gtex = l[1]
		rsid = l[0]
		rsid_snp[rsid] = gtex
		gtex_fst[gtex] = {}
		gtex_fst[gtex]['score'] = None
		gtex_fst[gtex]['pvalue'] = None
		gtex_fst[gtex]['rsid'] = None

with gzip.open(fst, 'rb') as f:
	next(f)
	for line in f:
		l = line.strip().split()
		rsid = l[0]
		score = float(l[3])
		pvalue = float(l[4])
		if rsid in rsid_snp:
			gtex_hit = rsid_snp[rsid]
			gtex_fst[gtex_hit]['score'] = score
			gtex_fst[gtex_hit]['pvalue'] = pvalue 
			gtex_fst[gtex_hit]['rsid'] = rsid

outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/gtex_snp_to_fst_map.tsv'
df = pandas.DataFrame.from_dict(gtex_fst, orient='index')
df.to_csv(outfile, sep='\t', index=True, header=True, index_label='variant_id')
