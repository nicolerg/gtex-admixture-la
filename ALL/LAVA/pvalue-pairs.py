# python script to get p-values from eGenes where lead SNP differs

import sys
import gzip 

subset = sys.argv[1]
allpairs = sys.argv[2]
outfile = sys.argv[3]

d = {}

with open(subset, 'rb') as sub:
	next(sub)
	for line in sub:
		l = line.strip().split()
		d[l[0]] = [ l[5], l[11] ]

header=True
with gzip.open(allpairs, 'rb') as allp, open(outfile, 'wb') as out:
	for line in allp:
		if header == True:
			out.write(line.strip()+'\tmethod_where_lead\n')
			header = False
			continue 
		l = line.strip().split()
		if l[0] in d:
			if l[1] == d[l[0]][0]:
				out.write(line.strip()+'\tglobal\n')
			elif l[1] == d[l[0]][1]:
				out.write(line.strip()+'\tlava\n')
