import gzip 
import pandas as pd
import sys 

geudir = "/mnt/lab_data/montgomery/nicolerg/geuvadis/"

afr_expr = geudir + "expression/GV-AFR.expression.bed.gz"
eur_expr = geudir + "expression/GV-EUR.expression.bed.gz"

# gene_snp_map = eqtldir + "/eur.eqtl.afr.0slope.txt"
# afr_freq = vcfdir + "/tg.AFR.0slope.freq"
# eur_vcf = vcfdir + "/tg.EUR.recode.0slope.recode.vcf.gz"

gene_snp_map = sys.argv[1]
freq = sys.argv[2]
vcf = sys.argv[3]
mono = sys.argv[4]

# make gene-snp map
gs = {}

with open(gene_snp_map, 'rb') as m:
	next(m)
	for line in m:
		l = line.strip().split()
		gs[l[0]] = l[1]

# get 0slope genotype at every snp
s = {}

with open(freq,'rb') as af:
	next(af)
	for line in af:
		l = line.strip().split()
		rsid = l[6]
		ref = l[4]
		alt = l[5]
		if '1' in ref:
			a = ref.split(':')[0]
		elif '1' in alt:
			a = alt.split(':')[0]
		geno = a+a 

		s[rsid] = geno

# read in expression
afrexpr = pd.read_csv(afr_expr, sep='\t', header=0)
eurexpr = pd.read_csv(eur_expr, sep='\t', header=0)

# read in poly genotypes
polyvcf = pd.read_csv(vcf, sep='\t', header=0)

eur_subjs = eurexpr.columns.values.tolist()[4:]
afr_subjs = afrexpr.columns.values.tolist()[4:]

# for every SNP-gene pair in the map, get EUR and AFR expression and genotypes
dflist = [None]*len(gs)
if mono == 'AFR':
	mono_subjs = afr_subjs
	poly_subjs = eur_subjs
	outfile = '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR-into-AFR/EUR.eqtl.0AFR.var.expr.tsv.gz'
	pop = ['EUR']*len(eur_subjs) + ['AFR']*len(afr_subjs)
elif mono == 'EUR':
	mono_subjs = eur_subjs
	poly_subjs = afr_subjs
	outfile = '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR/AFR.eqtl.0EUR.var.expr.tsv.gz'
	pop = ['AFR']*len(afr_subjs) + ['EUR']*len(eur_subjs)
j = 0
for gene in gs:
	variant = gs[gene]
	print gene + ': ' + variant
	# get EUR expression
	eur_lvls = eurexpr.loc[eurexpr['gene_id'] == gene, eur_subjs].values[0].tolist()
	# get AFR expression
	afr_lvls = afrexpr.loc[eurexpr['gene_id'] == gene, afr_subjs].values[0].tolist()
	# get AFR genotype

	mono_geno = [s[variant]]*len(mono_subjs)
	poly_geno = [None]*len(poly_subjs)
	
	# get poly genotype
	
	ref = polyvcf.loc[polyvcf['ID'] == variant, 'REF'].values[0]
	alt = polyvcf.loc[polyvcf['ID'] == variant, 'ALT'].values[0]
	for i in range(len(poly_subjs)):
		subj = poly_subjs[i]
		haps = polyvcf.loc[polyvcf['ID'] == variant, subj].values[0] # 0/1 | 0/1
		if haps == '1|1':
			poly_geno[i] = ref+ref
		elif haps == '1|0' or haps == '0|1':
			poly_geno[i] = ref+alt
		elif haps == '0|0':
			poly_geno[i] = alt+alt
		else:
			print 'Invalid genotype: subject ' + subj + ', variant: ' + variant

	# make a data frame
	df = pd.DataFrame(
		{'gene': [gene]*(len(poly_subjs)+len(mono_subjs)),
		'variant': [variant]*(len(poly_subjs)+len(mono_subjs)),
		'subjid': poly_subjs + mono_subjs,
		'pop': pop,
		'genotype': poly_geno + mono_geno,
		'expr': eur_lvls + afr_lvls
		})
	dflist[j] = df 
	j += 1

merged = pd.concat(dflist)

# write out table
merged.to_csv(outfile, sep='\t', index=False, compression='gzip')
