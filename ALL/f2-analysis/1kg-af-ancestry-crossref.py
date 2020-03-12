import gzip
import os
import sys

tgvcf_dir=sys.argv[1]
gtexla_dir=sys.argv[2]
outdir=sys.argv[3]
chr=sys.argv[4]
print 'Chromosome '+chr

tgvcf=tgvcf_dir+'/f2-1kg-intersect-chr'+chr+'.vcf.gz'
gtexla=gtexla_dir+'/bedintersect.chr'+chr+'.out'

outf=outdir+'/1kgAF-gtexLA-chr'+chr+'.tsv'

popdict={}

with open(gtexla,'r') as la:   
	next(la)
	for line in la:
		line = line.strip().split()
		popdict[int(line[1])] = line[7] # {start:population}

total = 0
matches = 0
afr_match = 0
asn_match = 0 
eur_match = 0
afr_count = 0
asn_count = 0
eur_count = 0
multia = 0
notfound = set()
with gzip.open(tgvcf,'r') as tg, open(outf, 'w') as out:
	out.write("#CHROM\tPOS\tID\tLOCAL_ANC\tTG_EAS_AF\tTG_EUR_AF\tTG_AFR_AF\tMATCH\n")
	for line in tg:
		if line[0] == "#":
			continue
		line = line.strip().split()
		pos = int(line[1])
		id = line[2]
		if pos not in popdict:
			notfound.add(pos)
			continue
		local_anc = popdict[pos]
		# get AFs
		info = line[7].strip().split(';')
		eas_af = [ x for x in info if 'EAS_AF' in x ][0].strip().split('=')[1]
		
		if ',' in eas_af: # multiallelic variant
			multia += 1
			continue

		if local_anc == 'AFR':
			afr_count += 1
		elif local_anc == 'ASN':
			asn_count += 1
		elif local_anc == 'EUR':
			eur_count += 1

		eas_af = float(eas_af)
		afr_af = float([ x for x in info if 'AFR_AF' in x ][0].strip().split('=')[1])
		eur_af = float([ x for x in info if 'EUR_AF' in x ][0].strip().split('=')[1])

		if local_anc == "AFR" and afr_af > eur_af and afr_af > eas_af:
			afr_match += 1
			match = 1
			matches += 1
		elif local_anc == "ASN" and eas_af > eur_af and eas_af > afr_af:
			asn_match += 1
			match = 1
			matches += 1
		elif local_anc == "EUR" and eur_af > afr_af and eur_af > eas_af:
			eur_match += 1
			match = 1
			matches += 1
		else:
			match = 0
		total += 1
		out.write('chr'+chr+'\t'+str(pos)+'\t'+id+'\t'+local_anc+'\t'+str(eas_af)+'\t'+str(eur_af)+'\t'+str(afr_af)+'\t'+str(match)+'\n')

print 'Fraction of total matches: '+ str(round(float(matches)/float(total),3))
print 'Fraction of African matches: '+ str(round(float(afr_match)/float(afr_count),3))
print 'Fraction of Asian matches: '+ str(round(float(asn_match)/float(asn_count),3))
print 'Fraction of European matches: '+ str(round(float(eur_match)/float(eur_count),3))
print 'Number of multiallelic variants ignored: ' + str(multia)
print 'Variants not found in BED file: ' + str(len(notfound)) + ' out of ' + str(total)
print 
