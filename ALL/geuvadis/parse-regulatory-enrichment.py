# DEPRECATED - specific for GEUVADIS overlaps

import gzip
import pandas as pd 
import sys
import os.path

indir = "/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Lung/chr1/test/"
gffdb = indir+'reg-enrich/homo_sapiens.GRCh38.Lung.Regulatory_Build.regulatory_activity.20161111.gff.gz'

globaleqtl = sys.argv[1]
#indir+'global.gtex.cutoff.eqtl.genes.tsv'
lavaeqtl = sys.argv[2]
#indir+'lava.gtex.cutoff.eqtl.genes.tsv'

if len(sys.argv) > 3:
	geneset = sys.argv[3]
	suffix = sys.argv[4]
else:
	geneset = 'none'
	suffix = 'all.eqtl.genes'

gset = False
if os.path.isfile(geneset):
	gset = True
	# make set 
	genes = set()
	with open(geneset, 'rb') as g:
		for line in g:
			l = line.strip().split('.')[0]
			genes.add(l)

globald={}
locald={}

def get_snps(file,d,genes=False):
	if file.endswith(".gz"):
		with gzip.open(file, 'rb') as gtex:
			next(gtex)
			for line in gtex:
				l = line.strip().split()
				if len(l) < 6:
					print l
					continue
				variant = l[5]
				gene = l[0].split('.')[0]
				if genes:
					if not gene in genes:
						continue
				print gene
				pos = int(float(variant.split('_')[1]))
				tss_dist = int(float(l[3]))
				chr = variant.split('_')[0]

				d[gene] = {}
				d[gene]['variant'] = variant
				d[gene]['chr'] = chr
				d[gene]['pos'] = pos
				d[gene]['tssdist'] = tss_dist
				d[gene]['CTCF_binding_site'] = 'NA'
				d[gene]['CTCF_binding_site_activity'] = 'NA'
				d[gene]['enhancer'] = 'NA'
				d[gene]['enhancer_activity'] = 'NA'
				d[gene]['open_chromatin_region'] = 'NA'
				d[gene]['promoter'] = 'NA'
				d[gene]['promoter_activity'] = 'NA'
				d[gene]['promoter_flanking_region'] = 'NA'
				d[gene]['promoter_flanking_region_activity'] = 'NA'
				d[gene]['TF_binding_site'] = 'NA'
	else:
		with open(file, 'rb') as gtex:
			next(gtex)
			for line in gtex:
				l = line.strip().split()
				variant = l[5]
				gene = l[0].split('.')[0]
				if genes:
					if not gene in genes:
						continue
				pos = int(float(variant.split('_')[1]))
				tss_dist = int(float(l[3]))
				chr = variant.split('_')[0]

				d[gene] = {}
				d[gene]['variant'] = variant
				d[gene]['chr'] = chr
				d[gene]['pos'] = pos
				d[gene]['tssdist'] = tss_dist
				d[gene]['CTCF_binding_site'] = 'NA'
				d[gene]['CTCF_binding_site_activity'] = 'NA'
				d[gene]['enhancer'] = 'NA'
				d[gene]['enhancer_activity'] = 'NA'
				d[gene]['open_chromatin_region'] = 'NA'
				d[gene]['promoter'] = 'NA'
				d[gene]['promoter_activity'] = 'NA'
				d[gene]['promoter_flanking_region'] = 'NA'
				d[gene]['promoter_flanking_region_activity'] = 'NA'
				d[gene]['TF_binding_site'] = 'NA'

	return d 

if gset:
	globald = get_snps(globaleqtl,globald,genes)
	locald = get_snps(lavaeqtl,locald,genes)
else:
	globald = get_snps(globaleqtl,globald)
	locald = get_snps(lavaeqtl,locald)

def compare_reg(d,start,stop,chr,egelm,ensembl,activity):
	for gene in d:
		if d[gene]['chr'] != chr:
			continue
		if (start <= d[gene]['pos']) and (stop >= d[gene]['pos']):
			d[gene][regelm] = ensembl
			if regelm != 'TF_binding_site' and regelm != 'open_chromatin_region':
				d[gene][regelm+'_activity'] = activity
	return d

# iterate through gff to annotate SNPs
with gzip.open(gffdb,'rb') as gff:
	for line in gff:
		l = line.strip().split('\t')
		chr = 'chr'+l[0]

		regelm = l[2]
		start = int(float(l[3]))
		stop = l[4]
		stop = int(float(l[4]))
		info = l[8].split(';')
		if len(info) < 7:
			print info
			continue
		activity = info[0].split('=')[1]
		ensembl = info[6].split('=')[1]

		globald = compare_reg(globald,start,stop,chr,regelm,ensembl,activity)
		locald = compare_reg(locald,start,stop,chr,regelm,ensembl,activity)

# write dictionaries out to file 
globaldf = pd.DataFrame.from_dict(globald,orient='index')
print globaldf
globaldf.to_csv(indir+'/reg-enrich/global.regenrich.'+suffix+'.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')
localdf = pd.DataFrame.from_dict(locald,orient='index')
print localdf
localdf.to_csv(indir+'/reg-enrich/lava.regenrich.'+suffix+'.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')

# # make a dictionary of LAVA/global eGene-SNP pairs
# inpairs = "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/aa-only/lung-LAVA-global-merged-allpairs-chr1.tsv.gz"
# # gene_id variant_id tss_distance ma_samples ma_count global_pval_nominal global_slope global_slope_se lava_pval_nominal lava_slope lava_slope_se afr_pval_nominal afr_slope afr_slope_se

# with gzip.open(inpairs, 'rb') as gtex:
# 	next(gtex)
# 	for line in gtex:
# 		l = line.strip().split()
# 		variant = l[1]
# 		pos = int(float(variant.split('_')[1]))
# 		tss_dist = int(float(l[2]))
# 		gene = l[0]
		
# 		if l[5] != 'NA' and l[5] != 'ERROR':
# 			global_pval_nominal = float(l[5])
# 		else:
# 			global_pval_nominal = 1.0
		
# 		if l[8] != 'NA' and l[8] != 'ERROR':
# 			lava_pval_nominal = float(l[8])
# 		else:
# 			lava_pval_nominal = 1.0

# 		if global_pval_nominal < cutoff:
# 			if gene in globalqtl:
# 				if globalqtl[gene]['pval'] > global_pval_nominal:
# 					globalqtl[gene]['pval'] = global_pval_nominal
# 					globalqtl[gene]['variant'] = variant
# 					globalqtl[gene]['pos'] = pos
# 					globalqtl[gene]['tssdist'] = tss_dist
# 			else:
# 				globalqtl[gene] = {}
# 				globalqtl[gene]['pval'] = global_pval_nominal
# 				globalqtl[gene]['variant'] = variant
# 				globalqtl[gene]['pos'] = pos
# 				globalqtl[gene]['tssdist'] = tss_dist
# 				globalqtl[gene]['CTCF_binding_site'] = 'NA'
# 				globalqtl[gene]['CTCF_binding_site_activity'] = 'NA'
# 				globalqtl[gene]['enhancer'] = 'NA'
# 				globalqtl[gene]['enhancer_activity'] = 'NA'
# 				globalqtl[gene]['open_chromatin_region'] = 'NA'
# 				globalqtl[gene]['promoter'] = 'NA'
# 				globalqtl[gene]['promoter_activity'] = 'NA'
# 				globalqtl[gene]['promoter_flanking_region'] = 'NA'
# 				globalqtl[gene]['promoter_flanking_region_activity'] = 'NA'
# 				globalqtl[gene]['TF_binding_site'] = 'NA'

# 		if lava_pval_nominal < cutoff:
# 			if gene in localqtl:
# 				if localqtl[gene]['pval'] > global_pval_nominal:
# 					localqtl[gene]['pval'] = global_pval_nominal
# 					localqtl[gene]['variant'] = variant
# 					localqtl[gene]['pos'] = pos
# 					localqtl[gene]['tssdist'] = tss_dist
# 			else:
# 				localqtl[gene] = {}
# 				localqtl[gene]['pval'] = global_pval_nominal
# 				localqtl[gene]['variant'] = variant
# 				localqtl[gene]['pos'] = pos
# 				localqtl[gene]['tssdist'] = tss_dist
# 				localqtl[gene]['CTCF_binding_site'] = 'NA'
# 				localqtl[gene]['CTCF_binding_site_activity'] = 'NA'
# 				localqtl[gene]['enhancer'] = 'NA'
# 				localqtl[gene]['enhancer_activity'] = 'NA'
# 				localqtl[gene]['open_chromatin_region'] = 'NA'
# 				localqtl[gene]['promoter'] = 'NA'
# 				localqtl[gene]['promoter_activity'] = 'NA'
# 				localqtl[gene]['promoter_flanking_region'] = 'NA'
# 				localqtl[gene]['promoter_flanking_region_activity'] = 'NA'
# 				localqtl[gene]['TF_binding_site'] = 'NA'

# # iterate through gff to annotate SNPs
# with gzip.open(gffdb,'rb') as gff:
# 	for line in gff:
# 		l = line.strip().split('\t')
# 		chr = 'chr'+l[0]

# 		if chr == 'chr1':
# 			regelm = l[2]
# 			start = int(float(l[3]))
# 			stop = l[4]
# 			stop = int(float(l[4]))
# 			info = l[8].split(';')
# 			if len(info) < 7:
# 				print info
# 				continue
# 			activity = info[0].split('=')[1]
# 			ensembl = info[6].split('=')[1]

# 			# first, iterate through global
# 			for gene in globalqtl:
# 				if (start <= globalqtl[gene]['pos']) and (stop >= globalqtl[gene]['pos']):
# 					globalqtl[gene][regelm] = ensembl
# 					if regelm != 'TF_binding_site' and regelm != 'open_chromatin_region':
# 						globalqtl[gene][regelm+'_activity'] = activity

# 			# next, iterate through local
# 			for gene in localqtl:
# 				if (start <= localqtl[gene]['pos']) and (stop >= localqtl[gene]['pos']):
# 					localqtl[gene][regelm] = ensembl
# 					if regelm != 'TF_binding_site' and regelm != 'open_chromatin_region':
# 						localqtl[gene][regelm+'_activity'] = activity

# # write dictionaries out to file 
# globaldf = pd.DataFrame.from_dict(globalqtl,orient='index')
# globaldf.to_csv(indir+'global-egene-regulatory-enrichment.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')
# localdf = pd.DataFrame.from_dict(localqtl,orient='index')
# localdf.to_csv(indir+'LAVA-egene-regulatory-enrichment.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')
