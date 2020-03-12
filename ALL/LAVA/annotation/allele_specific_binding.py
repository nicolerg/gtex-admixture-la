# parse allele-specific binding 

import sys 
import gzip 
import pandas
import cPickle as pickle
import os.path

annotation_dir = sys.argv[1]

asb = annotation_dir + '/annotation/allele-specific-binding/ng.3432-S7.gz'
rsid_map = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.rsid.map.txt.gz'

rsid_snp = {}
gtex_asb = {}

with gzip.open(rsid_map, 'rb') as snpmap:
	for line in snpmap:
		l = line.strip().split('\t')
		gtex = l[2]
		rsid = l[0]
		rsid_snp[rsid] = gtex
		gtex_asb[gtex] = 'NA'

with gzip.open(asb, 'rb') as asbs:
	next(asbs)
	for line in asbs:
		l = line.strip().split()
		rsids = l[3].strip().split(',')
		for x in rsids:
			if x in rsid_snp:
				gtex_hit = rsid_snp[x]
				gtex_asb[gtex_hit] = x

outfile = annotation_dir + '/allele-specific-binding/gtex_mac10_snp_to_asb_map.tsv'
df = pandas.DataFrame.from_dict(gtex_asb, orient='index')
df.to_csv(outfile, sep='\t', index=True, header=False)
pickle.dump(df, open(annotation_dir+'/allele-specific-binding/gtex_mac10_snp_to_asb_map.pkl', 'wb'))

