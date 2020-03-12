import gzip
import pandas as pd 
import sys

allp = sys.argv[1]
outfile = sys.argv[2]
indir = '/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/'

# make a dictionary of counts

gcounts = {}

with gzip.open(indir+allp, 'rb') as all:
	next(all)
	for line in all:
		l = line.strip().split()
		if l[0] in gcounts:
			gcounts[l[0]]['count'] += 1
			if float(l[3]) < gcounts[l[0]]['pval_nominal']:
				gcounts[l[0]]['variant_id'] = l[1]
				gcounts[l[0]]['pval_nominal'] = float(l[3])
				gcounts[l[0]]['slope'] = float(l[4])
				gcounts[l[0]]['tss_distance'] = int(l[2])
		else:
			gcounts[l[0]] = {}
			gcounts[l[0]]['count'] = 1
			gcounts[l[0]]['variant_id'] = l[1]
			gcounts[l[0]]['pval_nominal'] = float(l[3])
			gcounts[l[0]]['slope'] = float(l[4])
			gcounts[l[0]]['tss_distance'] = int(l[2])

egenes = pd.DataFrame.from_dict(gcounts,orient='index')
egenes.to_csv(indir+outfile,sep='\t',index=True,index_label='gene_id',compression='gzip')
