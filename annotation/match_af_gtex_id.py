import gzip 
import sys

outdir=sys.argv[1]

# convert BED output to GTEx variant_id format (called by get_1kg_af.sh)

info = outdir + '/gtex.admixed.merged.1kg.AF.INFO'
hg38 = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/gtex.admixed.MAC10.all.snps.txt.gz'

# CHROM   POS     REF     ALT     EAS_AF  EUR_AF  AFR_AF
# 1       858792  G       A,T     0,0     0.007,0 0.0008,0.1967
# 1       58814   G       A       0.0228  0.0676  0.1967

snpmap = {} # snpmap[chr1][21243] = chr1_21243_A_C_b38
for c in range(1,23):
	snpmap['chr'+str(c)] = {}

with gzip.open(hg38, 'rb') as b38:
	for line in b38:
		l = line.strip().split('_')
		chrom = l[0]
		pos = int(l[1])
		ref = l[2]
		alt = l[3]
		if chrom in snpmap:
			snpmap[chrom][pos] = line.strip()

outfile = outdir+ '/gtex.admixed.MAC10.snps.af.txt'
# CHROM   POS     REF     ALT     EAS_AF  EUR_AF  AFR_AF
with open(info, 'rb') as af, open(outfile, 'wb') as out:
	out.write('\t'.join(['CHROM','POS','REF','ALT','EAS_AF','EUR_AF','AFR_AF','variant_id'])+'\n')
	next(af)
	for line in af:
		l = line.strip().split()
		chrom = 'chr'+l[0]
		
		if chrom in snpmap:
			pos = int(l[1])
			
			if pos in snpmap[chrom]:
				variant_id = snpmap[chrom][pos]
				gtex_ref = variant_id.split('_')[2]
				gtex_alt = variant_id.split('_')[3]
				kg_ref = l[2]
				kg_alt = l[3]

				if "," in kg_alt:
					kg_alt_list = kg_alt.split(',')
					for i in range(len(kg_alt_list)):
						alt = kg_alt_list[i]
						if gtex_ref == kg_ref and gtex_alt == alt:
							eas_af = float(l[4].split(',')[i])
							eur_af = float(l[5].split(',')[i])
							afr_af = float(l[6].split(',')[i])
							out.write('\t'.join([chrom, str(pos), gtex_ref, gtex_alt, str(eas_af), str(eur_af), str(afr_af), variant_id]) + '\n')
							break

				if gtex_ref == kg_ref and gtex_alt == kg_alt:
					out.write(line.strip() + '\t' + variant_id + '\n')
				
				elif gtex_ref == kg_alt and gtex_alt == kg_ref:
					eas_af = 1 - float(l[4])
					eur_af = 1 - float(l[5])
					afr_af = 1 - float(l[6])
					out.write('\t'.join([chrom, str(pos), gtex_ref, gtex_alt, str(eas_af), str(eur_af), str(afr_af), variant_id]) + '\n')
