# CADD 

import sys 
import gzip 
import pandas
import pickle

# ## CADD GRCh38-v1.4 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2018. All
# #Chrom  Pos     Ref     Alt     RawScore        PHRED
# 1       10001   T       A       0.591814        8.493

all_snps = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.all.snps.txt.gz' 
cadd = '/mnt/lab_data/montgomery/nicolerg/RESOURCES/cadd/whole_genome_SNVs.tsv.gz' #hg38 

snpmap = {}
for c in range(1,23):
	snpmap['chr'+str(c)] = {}

master = {} # gtex_snp:CADD_score

with gzip.open(all_snps, 'rb') as gtex:
	for line in gtex:
		l = line.strip().split('_')
		chrom = l[0]
		if chrom in snpmap:
			snpmap[chrom][int(l[1])] = {}
			snpmap[chrom][int(l[1])]['snp'] = line.strip()
			snpmap[chrom][int(l[1])]['ref'] = l[2]
			snpmap[chrom][int(l[1])]['alt'] = l[3]

with gzip.open(cadd, 'rb') as cad:
	for line in cad:
		if line.startswith('#'):
			continue
		l = line.strip().split()
		chrom = 'chr'+l[0]
		loc = int(l[1])
		ref = l[2]
		alt = l[3]
		score = float(l[4])
		if chrom in snpmap:
			if loc in snpmap[chrom]:
				if snpmap[chrom][loc]['ref'] == ref and snpmap[chrom][loc]['alt'] == alt:
					name = snpmap[chrom][loc]['snp']
					master[name] = score
				elif snpmap[chrom][loc]['ref'] == alt and snpmap[chrom][loc]['alt'] == ref:
					name = snpmap[chrom][loc]['snp']
					master[name] = score

# write a map of CADD scores for all relevant SNPs

with open('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cadd/gtex.admixed.MAC10.all.snps_cadd.pkl', 'wb') as handle:
    pickle.dump(master, handle, protocol=pickle.HIGHEST_PROTOCOL)

outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cadd/gtex.admixed.MAC10.all.snps_cadd.tsv'
df = pandas.DataFrame.from_dict(master, orient='index')
df.to_csv(outfile, sep='\t', index=True, header=False)
