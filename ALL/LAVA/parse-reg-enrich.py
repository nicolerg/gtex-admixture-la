import gzip
import pandas as pd 
import sys
import os.path

merged = sys.argv[1]
gffdb = sys.argv[2]
outdir = sys.argv[3]

global_eqtl={}
lava_eqtl={}
global_snp={}
lava_snp={}

def fill_dict(d, gene, variant, tss_dist):
	pos = int(float(variant.split('_')[1]))
	chrom = variant.split('_')[0]
	d[gene] = {}
	d[gene]['variant'] = variant
	d[gene]['chr'] = chrom
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

def get_eSNP(line,global_snp,lava_snp):
	l = line.strip().split()
	global_variant = l[5]
	lava_variant = l[11]
	gene = l[0].split('.')[0]
	global_snp = fill_dict(global_snp, gene, global_variant, 0)
	lava_snp = fill_dict(lava_snp, gene, lava_variant, 0)
	return [global_snp,lava_snp]

def get_snp(line,global_eqtl,lava_eqtl,global_snp,lava_snp):
	l = line.strip().split()
	global_variant = l[5]
	lava_variant = l[10]
	gene = l[0].split('.')[0]
	global_dist = int(float(l[4]))
	lava_dist = int(float(l[9]))
	# keep track of only eQTLs
	if float(l[6]) < float(l[1]):
		global_eqtl = fill_dict(global_eqtl, gene, global_variant, global_dist)
	if float(l[11]) < float(l[1]):
		lava_eqtl = fill_dict(lava_eqtl, gene, lava_variant, lava_dist)
	# keep track of all lead SNPs
	global_snp = fill_dict(global_snp, gene, global_variant, global_dist)
	lava_snp = fill_dict(lava_snp, gene, lava_variant, lava_dist)
	return [global_eqtl,lava_eqtl,global_snp,lava_snp]

def compare_reg(d,start,stop,chrom,regelm,ensembl,activity):
	for gene in d:
		if d[gene]['chr'] != chrom:
			continue
		if (start <= d[gene]['pos']) and (stop >= d[gene]['pos']):
			d[gene][regelm] = ensembl
			if regelm != 'TF_binding_site' and regelm != 'open_chromatin_region':
				d[gene][regelm+'_activity'] = activity
	return d

# with open(merged,'rb') as egenes:
# 	next(egenes)
# 	for line in egenes:
# 		global_eqtl,lava_eqtl,global_snp,lava_snp = get_snp(line,global_eqtl,lava_eqtl,global_snp,lava_snp)

with open(merged,'rb') as egenes:
	next(egenes)
	for line in egenes:
		global_snp,lava_snp = get_eSNP(line,global_snp,lava_snp)

# iterate through gff to annotate SNPs
with gzip.open(gffdb,'rb') as gff:
	for line in gff:
		l = line.strip().split('\t')
		chrom = 'chr'+l[0]

		regelm = l[2]
		start = int(float(l[3]))
		stop = int(float(l[4]))
		info = l[8].split(';')
		if len(info) < 7:
			print info
			continue
		activity = info[0].split('=')[1]
		ensembl = info[6].split('=')[1]

		#global_eqtl = compare_reg(global_eqtl,start,stop,chrom,regelm,ensembl,activity)
		#lava_eqtl = compare_reg(lava_eqtl,start,stop,chrom,regelm,ensembl,activity)
		global_snp = compare_reg(global_snp,start,stop,chrom,regelm,ensembl,activity)
		lava_snp = compare_reg(lava_snp,start,stop,chrom,regelm,ensembl,activity)

# # write dictionaries out to file 
# global_eqtl_df = pd.DataFrame.from_dict(global_eqtl,orient='index')
# global_eqtl_df.to_csv(outdir+'/global.eqtl.regenrich.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')

# global_snp_df = pd.DataFrame.from_dict(global_snp,orient='index')
# global_snp_df.to_csv(outdir+'/global.leadsnp.regenrich.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')

# lava_eqtl_df = pd.DataFrame.from_dict(lava_eqtl,orient='index')
# lava_eqtl_df.to_csv(outdir+'/lava.eqtl.regenrich.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')

# lava_snp_df = pd.DataFrame.from_dict(lava_snp,orient='index')
# lava_snp_df.to_csv(outdir+'/lava.leadsnp.regenrich.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')

# write dictionaries out to file 

global_snp_df = pd.DataFrame.from_dict(global_snp,orient='index')
global_snp_df.to_csv(outdir+'/global.merged.eQTL.diff.leadSNP.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')

lava_snp_df = pd.DataFrame.from_dict(lava_snp,orient='index')
lava_snp_df.to_csv(outdir+'/lava.merged.eQTL.diff.leadSNP.tsv.gz',sep='\t',index=True,index_label='gene',compression='gzip')
