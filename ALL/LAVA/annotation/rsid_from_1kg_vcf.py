#!/bin/python

import gzip
import sys 
import pandas 

# make RSID map from 1KG. this is not comprehensive, but it's what we need to get Fst values 

chrom = sys.argv[1]

infile = '/mnt/lab_data/montgomery/shared/1KG/hg38/ALL.chr{0}_GRCh38.genotypes.20170504.vcf.gz'.format(chrom)
all_snps = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.all.snps.txt.gz'

rsid_to_snp = {} # rsid:gtex_snp (1 to 1)

chrom = 'chr' + chrom

snpmap = {}
rsid_to_snp = {}

with gzip.open(all_snps, 'rb') as gtex:
	for line in gtex:
		l = line.strip().split('_')
		c = l[0]
		if c == chrom:
			position = int(l[1])
			alt = l[3]
			if position not in snpmap:
				snpmap[position] = {}
			snpmap[position][alt] = line.strip() # chr10_45349_G_A_b38 --> 45349: { A: chr10_45349_G_A_b38 , T: chr10_45349_G_T_b38 } 

with gzip.open(infile, 'rb') as vcf:
	for line in vcf:
		if line.startswith('#'):
			continue
		l = line.strip().split()
		if line.startswith('chr'):
			c = l[0]
		else:
			c = 'chr'+l[0]
		if c == chrom:

			pos = int(l[1])

			if pos in snpmap:

				rsid = l[2]
				kg_ref = l[3]
				kg_alt = l[4]

				if kg_alt in snpmap[pos]:

					gtex_id = snpmap[pos][kg_alt]
					gtex_other = gtex_id.split('_')[2]
					if gtex_other == kg_ref:
						rsid_to_snp[rsid] = gtex_id
						print gtex_id

				elif kg_ref in snpmap[pos]: 

					gtex_id = snpmap[pos][kg_ref]
					gtex_other = gtex_id.split('_')[3]
					if gtex_other == kg_alt:
						rsid_to_snp[rsid] = gtex_id
						print gtex_id

outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.rsid_map_from_1kg_{}.txt.gz'.format(chrom)
df = pandas.DataFrame.from_dict(rsid_to_snp, orient='index')
df.to_csv(outfile, sep='\t', index=True, index_label='rsid', header=True, compression='gzip')
