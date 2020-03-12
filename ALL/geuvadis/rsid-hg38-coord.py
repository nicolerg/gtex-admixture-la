import gzip
import cPickle as pickle
import os

indir = '/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/'
outbed = '/users/nicolerg/gtex-admix/geuvadis.norsid.bed'

# if not os.path.isfile(indir+'geu.id.pkl') or not os.path.isfile(outbed) or not os.path.isfile(indir+'geuvadis.hg38.ids.bed'):

# 	# make a dictionary of all of the geuvadis rsids; write non-rsids to BED file
# 	geuid = {}
# 	with gzip.open(indir+'geuvadis.rsids.gz','rb') as rsid, open(outbed, 'wb') as bed:
# 		for line in rsid:
# 			line = line.strip()
# 			if line.startswith("rs"):
# 				geuid[line] = ''
# 			elif line.startswith("snp") or line.startswith("indel"):
# 				l = line.split('_')
# 				chr = "chr"+l[1].strip().split('_')[0]
# 				start = int(l[2]) - 1
# 				stop = int(l[2])
# 				bed.write('\t'.join([chr,str(start),str(stop),line])+'\n')
# 			elif line.startswith("esv"):
# 				continue

# 	with open(indir+'geu.id.pkl', 'wb') as handle:
# 	    pickle.dump(geuid, handle, protocol=2)

# 	exit()

# 	# outside of script, liftOver bed file --> geuvadis.hg38.ids.bed

## Previous steps:
# Extract all GEUVADIS IDs
# For IDs with a valid SNP pos (not ending in .5), write bed file
# For IDs with no valud SNP pos (~10k), cross-reference with dbsnp to get hg19 coordinates (parse-snp-hg19.py)
# Concatenate bed files
# liftOver

# with liftOver bed file, make dictionary of IDs and hg38 coordinates
if not os.path.isfile(indir+'geu.id.hg38.pkl'):

	geuid = {}
	with gzip.open(indir+'geuvadis.variants.hg38.bed.gz','rb') as rsid:
		for line in rsid:
			l = line.strip().split()
			chr = l[0]
			pos = l[2]
			id = l[3]
			geuid[id] = chr+':'+pos

	with open(indir+'geu.id.hg38.pkl', 'wb') as handle:
	    pickle.dump(geuid, handle, protocol=2)

# if not os.path.isfile(indir+'geu.hg38.id.pkl'):

# 	# add lifted over variants to dictionary

# 	with open(indir+'geu.id.pkl', 'rb') as handle:
# 	    geuid = pickle.load(handle)

# 	with open(indir+'geuvadis.hg38.ids.bed', 'rb') as bed:
# 		for line in bed:
# 			l = line.strip().split()
# 			chr = l[0].strip().split('_')[0]
# 			pos = l[2]
# 			id = l[3]
# 			geuid[id] = chr+':'+pos

# 	# iterate through snpdb and add hg38 coordinates
# 	with gzip.open(indir+'snp150All.bed.gz', 'rb') as snpdb:
# 		for line in snpdb:
# 			l = line.strip().split()
# 			if l[3] in geuid and not '_' in l[0]:
# 				geuid[l[3]] = l[0]+':'+l[2]

# 	with open(indir+'geu.hg38.id.pkl', 'wb') as handle:
# 	    pickle.dump(geuid, handle, protocol=2)

if not os.path.isfile(indir+'YRI89.gene.cis.FDR5.best.rs137.hg38.txt.gz') and os.path.isfile(indir+'geu.id.hg38.pkl'):

	with open(indir+'geu.id.hg38.pkl', 'rb') as handle:
	    geuid = pickle.load(handle)

	# # iterate through dictionary and remove keys with no hg38 id
	# count = 0
	# for k, v in geuid.items():
	# 	if v == '':
	# 		count += 1
	# 		print k
	# 		del geuid[k]

	# print count
	# print len(geuid)

	# exit()

	# iterate through geuvadis files and add hg38 ids
	header = 'SNP_ID\tID\tGENE_ID\tPROBE_ID\tCHR_SNP\tCHR_GENE\tSNPpos\tTSSpos\tdistance\trvalue\tpvalue\tlog10pvalue\thg38_id\n'
	geufiles = ['EUR373.gene.cis.FDR5.all.rs137.txt.gz','EUR373.gene.cis.FDR5.best.rs137.txt.gz','YRI89.gene.cis.FDR5.all.rs137.txt.gz','YRI89.gene.cis.FDR5.best.rs137.txt.gz']
	for f in geufiles:
		print f
		missing = 0
		outfile = '.'.join(f.split('.')[0:6])+'.hg38.txt.gz'
		with gzip.open(f,'rb') as geu, gzip.open(outfile, 'wb') as hg38:
			hg38.write(header)
			for line in geu:
				if line.startswith("SNP_ID"):
					continue
				l = line.strip().split()
				if l[0] in geuid:
					hg38_id = geuid[l[0]]
					hg38.write(line.strip() + '\t' + hg38_id + '\n')
				else:
					missing += 1
		print str(missing) + ' variants not lifted over to hg38'

if os.path.isfile(indir+'EUR373.gene.cis.FDR5.best.rs137.hg38.txt.gz') and not os.path.isfile(indir+'geu.eur.eqtl.pkl'):

	# for EUR eQTLs from GEU, what is the replication in GTEx tissues?
	# % of SNP-gene pairs with pval < 0.05, ordered by tissue --> write out to table for plotting in R
	# repeat for African eQTLs

	# make a dictionary of GEU EUR eQTLs (1 per gene)

	eurqtl = {}
	with gzip.open(indir+'EUR373.gene.cis.FDR5.best.rs137.hg38.txt.gz') as eur:
		next(eur)
		for line in eur:
			l = line.strip().split()
			hg38_id = l[12]
			genebase = l[2].strip().split('.')[0]
			eurqtl[hg38_id] = genebase
	print 'Number of EUR eQTLs after processing: ' + str(len(eurqtl))

	# make a dictionary of GEU AFR eQTLs (1 per gene)

	afrqtl = {}
	with gzip.open(indir+'YRI89.gene.cis.FDR5.best.rs137.hg38.txt.gz') as afr:
		next(afr)
		for line in afr:
			l = line.strip().split()
			hg38_id = l[12]
			genebase = l[2].strip().split('.')[0]
			afrqtl[hg38_id] = genebase
	print 'Number of AFR eQTLs after processing: ' + str(len(afrqtl))

	# make a dictionary of sample size
	
	samp_size = {}
	exprdir = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/'
	expr = os.listdir(exprdir)
	for f in expr:
		if not f.endswith('tbi'):
			tissue = f.strip().split('.')[0]
			with gzip.open(exprdir+f) as inexpr:
				l = inexpr.readline()
				nsamp = len(l.strip().split()) - 4
				samp_size[tissue] = nsamp

	with open(indir+'geu.eur.eqtl.pkl', 'wb') as handle:
		pickle.dump(eurqtl, handle, protocol=2)

	with open(indir+'geu.afr.eqtl.pkl', 'wb') as handle:
		pickle.dump(afrqtl, handle, protocol=2)

	with open(indir+'gtex.sampsize.pkl', 'wb') as handle:
		pickle.dump(samp_size, handle, protocol=2)


	# # for each GTEx tissue, count how many replicate with pval < 0.05

	# gtexdir = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/'

	# rep = {}

	# allpairs = os.listdir(gtexdir)
	# for f in allpairs:
	# 	afrcount = 0
	# 	eurcount = 0
	# 	tissue = f.strip().split('.')[0]
	# 	print tissue

	# 	with gzip.open(gtexdir+f, 'rb') as gtex:
	# 		next(gtex)
	# 		for line in gtex:
	# 			l = line.strip().split()
	# 			pval = float(l[6])
	# 			if pval < 0.05:
	# 				genebase = l[0].split('.')[0]
	# 				hg38_id = ':'.join(l[1].split('_')[0:2])
	# 				if hg38_id in eurqtl:
	# 					if eurqtl[hg38_id] == genebase:
	# 						eurcount += 1
	# 				if hg38_id in afrqtl:
	# 					if afrqtl[hg38_id] == genebase:
	# 						afrcount += 1

	# 	eurfrac = float(eurcount/len(eurqtl))
	# 	afrfrac = float(afrcount/len(afrqtl))

	# 	print "Eur rep: " + str(eurfrac)
	# 	print "Afr rep: " + str(afrfrac)

	# 	rep[tissue] = {}
	# 	rep[tissue]['EUR_frac'] = eurfrac
	# 	rep[tissue]['AFR_frac'] = afrfrac
	# 	rep[tissue]['sample_size'] = samp_size[tissue]

	# 	if counter == 3:
	# 		break

	# # convert dictionary to tsv; write to file
	# df = pd.DataFrame.from_dict(rep, orient='index')
	# df.to_csv(outdir+'geu.eur.eqtl.gtex.replication.bytissue.tsv', sep='\t',header=True,index=True)
