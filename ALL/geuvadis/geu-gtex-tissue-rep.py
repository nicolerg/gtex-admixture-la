import cPickle as pickle
import gzip
import os
import sys
import time

infile = sys.argv[1]
outdir = sys.argv[2] + '/'
indir = outdir
tissue = infile.strip().split('.')[0]
gtexdir = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/'

with open(indir+'geu.eur.eqtl.pkl', 'rb') as handle:
	eurqtl = pickle.load(handle)

with open(indir+'geu.afr.eqtl.pkl', 'rb') as handle:
	afrqtl = pickle.load(handle)

with open(indir+'gtex.sampsize.pkl', 'rb') as handle:
	samp_size = pickle.load(handle)

# count how many eqtl genes replicate with pval < 0.05

eurcount = 0
afrcount = 0

with gzip.open(gtexdir+infile, 'rb') as gtex:
	next(gtex)
	for line in gtex:
		l = line.strip().split()
		pval = float(l[6])
		if pval < 0.05:
			genebase = l[0].split('.')[0]
			hg38_id = ':'.join(l[1].split('_')[0:2])
			if hg38_id in eurqtl:
				if eurqtl[hg38_id] == genebase:
					eurcount += 1
					print hg38_id
			if hg38_id in afrqtl:
				if afrqtl[hg38_id] == genebase:
					afrcount += 1
					print hg38_id

eurfrac = float(eurcount)/float(len(eurqtl))
afrfrac = float(afrcount)/float(len(afrqtl))

print "Eur rep: " + str(eurfrac)
print "Afr rep: " + str(afrfrac)

n = samp_size[tissue]

with open(outdir + 'geu.gtex.replication.' + tissue + '.txt', 'wb') as out:
	#out.write('\t'.join(['tissue','EUR_frac','AFR_frac','sample_size']) + '\n')
	out.write('\t'.join([tissue, str(eurfrac), str(afrfrac), str(n)]) + '\n')
