import pandas as pd 
import sys

infile = sys.argv[1]
outdir = sys.argv[2]

altalt=[]
refref=[]
het=[]

with open(infile,'rb') as tstat:
	next(tstat)
	for line in tstat:
		l = line.strip().split()
		mono = l[3]
		poly = l[2]
		tstat = float(l[4])
		if len(mono) == 2:
			ref = mono[0]
		else:
			ref = mono[0:len(mono)/2]
		if poly == ref+ref:
			refref.append(tstat)
		elif poly.endswith(ref) or poly.startswith(ref):
			het.append(tstat)
		else:
			altalt.append(tstat)

label = ['altalt']*len(altalt) + ['refref']*len(refref) + ['het']*len(het)
values = altalt + refref + het

tstat_out = pd.DataFrame({'label': label,'tstat': values})
tstat_out.to_csv(outdir+'/tstat-per-genotype.tsv', sep='\t', header=True, index=False)

