# make GTEx SNP: RSID map

import sys 
import gzip 
import pandas

db38 = '/mnt/lab_data/montgomery/nicolerg/RESOURCES/dbsnp/snp150.hg38.txt.gz'
db19 = '/mnt/lab_data/montgomery/nicolerg/RESOURCES/dbsnp/snp150.hg19.txt.gz'
all_snps = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.all.snps.txt.gz'

rsid_to_snp = {} # rsid:gtex_snp (1 to 1)

snpmap = {}
for c in range(1,23):
	snpmap['chr'+str(c)] = {}

with gzip.open(all_snps, 'rb') as gtex:
	for line in gtex:
		l = line.strip().split('_')
		chrom = l[0]
		position = int(l[1])
		alt = l[3]
		if chrom in snpmap:
			if position not in snpmap[chrom]:
				snpmap[chrom][position] = {}
			snpmap[chrom][position][alt] = line.strip() # chr10_45349_G_A_b38 --> chr10: { 45349: { A: chr10_45349_G_A_b38 , T: chr10_45349_G_T_b38 } }

with gzip.open(db38, 'rb') as dbsnp:
	for line in dbsnp:
		l = line.strip().split()
		chrom = l[1]
		loc = int(l[3])
		rsid = l[4]

		alleles = l[9].split('/')

		if chrom not in snpmap:
			continue
		if loc in snpmap[chrom]:
			
			# check if it's the right alt/ref allele 

			for gtex_alt in snpmap[chrom][loc]:

				print gtex_alt
				print alleles

				if gtex_alt in alleles:

					rsid_to_snp[rsid] = {}
					name = snpmap[chrom][loc][gtex_alt]
					rsid_to_snp[rsid]['hg38'] = name
					break


# make a map from RSID to snp

with gzip.open(db19, 'rb') as b19:
	for line in b19:
		l = line.strip().split()
		chrom = l[1]
		loc = int(l[3])
		rsid = l[4]
		if rsid in rsid_to_snp:
			rsid_to_snp[rsid]['hg19'] = chrom+'_'+str(loc)

outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.rsid.map_fixed.txt.gz'
df = pandas.DataFrame.from_dict(rsid_to_snp, orient='index')
df.to_csv(outfile, sep='\t', index=True, index_label='rsid', header=True, compression='gzip')
