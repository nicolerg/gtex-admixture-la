import gzip
import pandas as pd 

outdir = "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/"
ids = "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/rsid.nopos.uniq"
snp19 = "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/snp150.hg19.txt.gz"

# make a dictionary of ids that need to be given positions

rsbed = {}
with open(ids, 'rb') as rsid:
	for line in rsid:
		rsid = line.strip()
		rsbed[rsid] = {}

# iterate through database
with gzip.open(snp19, 'rb') as snp:
	for line in snp:
		l = line.strip().split()
		if l[4] in rsbed:
			rsbed[l[4]]['CHR'] = l[1].split('_')[0]
			rsbed[l[4]]['START'] = l[2]
			rsbed[l[4]]['STOP'] = l[3]

# convert dictionary to tsv; write to file
df = pd.DataFrame.from_dict(rsbed, orient='index')
df.to_csv(outdir+'geu.rsid.nopos.hg19.tsv', sep='\t',header=True,index=True)
