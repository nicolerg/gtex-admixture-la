import gzip 
import pandas as pd
import sys 

geudir = "/mnt/lab_data/montgomery/nicolerg/geuvadis/"

afr_expr = geudir + "expression/GV-AFR.expression.bed.gz"
eur_expr = geudir + "expression/GV-EUR.expression.bed.gz"

gene_snp_map = "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/overlapping/overlapping.eqtls"
afr_vcf = "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/overlapping/tg.AFR.recode.overlap.recode.vcf.gz"
eur_vcf = "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/overlapping/tg.EUR.recode.overlap.recode.vcf.gz"

# make gene-snp map
gs = {}

with open(gene_snp_map, 'rb') as m:
	next(m)
	for line in m:
		l = line.strip().split()
		gs[l[0]] = l[1]

# read in expression
afrexpr = pd.read_csv(afr_expr, sep='\t', header=0)
eurexpr = pd.read_csv(eur_expr, sep='\t', header=0)

# read in genotypes
afr_vcf = pd.read_csv(afr_vcf, sep='\t', header=0)
eur_vcf = pd.read_csv(eur_vcf, sep='\t', header=0)

eur_subjs = eurexpr.columns.values.tolist()[4:]
afr_subjs = afrexpr.columns.values.tolist()[4:]

# for every SNP-gene pair in the map, get EUR and AFR expression and genotypes
outfile = "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/overlapping/overlapping.snp.expression.tsv.gz"

def get_geno(vcf,variant,subjects,pop,lvls):
	ref = vcf.loc[vcf['ID'] == variant, 'REF'].values[0]
	alt = vcf.loc[vcf['ID'] == variant, 'ALT'].values[0]
	geno = [None]*len(subjects)
	for i in range(len(subjects)):
		subj = subjects[i]
		haps = vcf.loc[vcf['ID'] == variant, subj].values[0] # 0/1 | 0/1
		if haps == '1|1':
			geno[i] = ref+ref
		elif haps == '1|0' or haps == '0|1':
			geno[i] = ref+alt
		elif haps == '0|0':
			geno[i] = alt+alt
		else:
			print 'Invalid genotype: subject ' + subj + ', variant: ' + variant

	# make a data frame
	df = pd.DataFrame(
		{'gene': [gene]*(len(subjects)),
		'variant': [variant]*(len(subjects)),
		'subjid': subjects,
		'pop': [pop]*(len(subjects)),
		'genotype': geno,
		'expr': lvls
		})
	return df

dflist = [None]*(len(gs)*2)
j = 0
for gene in gs:
	variant = gs[gene]
	print gene + ': ' + variant
	# get EUR expression
	eur_lvls = eurexpr.loc[eurexpr['gene_id'] == gene, eur_subjs].values[0].tolist()
	# get AFR expression
	afr_lvls = afrexpr.loc[eurexpr['gene_id'] == gene, afr_subjs].values[0].tolist()
	
	# get AFR genotype
	df = get_geno(afr_vcf,variant,afr_subjs,'AFR',afr_lvls)
	dflist[j] = df 
	j += 1
	# get EUR genotpe
	df = get_geno(eur_vcf,variant,eur_subjs,'EUR',eur_lvls)
	dflist[j] = df 
	j += 1

merged = pd.concat(dflist)

# write out table
merged.to_csv(outfile, sep='\t', index=False, compression='gzip')
