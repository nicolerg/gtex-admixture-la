"""
Run this script on the merged allpairs file output by combined lava & global eQTL
calling to summarize the top SNPs for each gene. All tied top eSNPs are output.

Currently, top SNPs for all tested genes are reported. They are not filtered by
significance. 

Updated 20 May 2019 
"""

import gzip
import pandas as pd 
import sys
import pickle 
import os 

def add_entry(d, line, pval_index, slope_index):
	l = line.strip().split()
	pval = float(l[pval_index])
	slope = float(l[slope_index])
	
	if l[0] in d:
		
		d[l[0]]['count'] += 1
		
		# break if pvalue is greater than current min 
		if pval > d[l[0]]['pval_nominal'][0]:
			return d

		# update dict if pvalue is less than current min
		if pval < d[l[0]]['pval_nominal']:
			d[l[0]]['variant_id'] = [ l[1] ]
			d[l[0]]['pval_nominal'] = [ pval ]
			d[l[0]]['slope'] = [ slope ]
			return d

		# append to existing entry if pvalues are the same
		if pval == d[l[0]]['pval_nominal']:
			d[l[0]]['variant_id'].append(l[1])
			d[l[0]]['slope'].append(slope)
			d[l[0]]['slope'].append(pval)
			return d

	else: # new entry
		
		d[l[0]] = {}
		d[l[0]]['count'] = 1
		d[l[0]]['variant_id'] = [ l[1] ]
		d[l[0]]['pval_nominal'] = [ pval ]
		d[l[0]]['slope'] = [ slope ]
		return d


def get_lead_snps(file, pval_index, slope_index):

	egenes = {}
	fail = [ 'NaN', 'ERROR', 'NA']

	with gzip.open(file, 'rb') as allpairs:
		next(allpairs)
		for line in allpairs:
			l = line.strip().split()
			
			if l[pval_index] in fail or l[slope_index] in fail:
				continue

			egenes = add_entry(egenes, line, pval_index, slope_index)

	return egenes

def dict_to_long_data(d):
	first = True
	for gene in d:

		df = pd.DataFrame.from_dict(d[gene])
		df['gene_id'] = gene

		if first:
			stretched = df
			first = False
		else:
			stretched = pd.concat([stretched, df])
	
	return stretched 

def parse_newpeers(tissue):

	allanc_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/both_anc/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-merged-allpairs-chr19.tsv.gz'
	newpeers = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/new_peers/eqtl/Skin_Not_Sun_Exposed_Suprapubic/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-allpairs-chr19.tsv.gz'

	# make egene dictionaries for each file 

	admixed_local_egenes = get_lead_snps(allanc_allpairs, 8, 9)
	admixed_global_egenes = get_lead_snps(allanc_allpairs, 5, 6)
	newpeer_global_egenes = get_lead_snps(newpeers, 5, 6)
	newpeer_local_egenes = get_lead_snps(newpeers, 8, 9)

	# save dictionaries 

	outdir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons'
	dicts = [ admixed_local_egenes, admixed_global_egenes, newpeer_local_egenes, newpeer_global_egenes ]
	names = [ '/both_anc/admixed_local_egenes', '/both_anc/admixed_global_egenes', '/new_peers/eqtl/Skin_Not_Sun_Exposed_Suprapubic/newpeer_local_egenes', '/new_peers/eqtl/Skin_Not_Sun_Exposed_Suprapubic/newpeer_global_egenes' ]

	for i in range(len(dicts)):
		with open(outdir + names[i] + '.' + tissue + '.pkl', 'wb') as handle:
			pickle.dump(dicts[i], handle, protocol=pickle.HIGHEST_PROTOCOL)

	return [admixed_local_egenes, admixed_global_egenes, newpeer_local_egenes, newpeer_global_egenes]



def parse_files(tissue):

	if tissue == 'Muscle_Skeletal':
		#my_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Muscle_Skeletal/Muscle_Skeletal-LAVA-global-merged-allpairs-chr19.tsv.gz'
		gtex_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/gtex_v8/Muscle_Skeletal.allpairs.txt.gz'
		allanc_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/both_anc/Muscle_Skeletal-LAVA-global-merged-allpairs-chr19.tsv.gz'
		allsubj_allpairs = False

	if tissue == 'Skin_Not_Sun_Exposed_Suprapubic':
		#my_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Skin_Not_Sun_Exposed_Suprapubic/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-merged-allpairs-chr19.tsv.gz'
		gtex_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/gtex_v8/Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz'
		allanc_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/both_anc/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-merged-allpairs-chr19.tsv.gz'
		allsubj_allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/all_ind/Skin_Not_Sun_Exposed_Suprapubic-LAVA-global-merged-allpairs-chr19.tsv.gz'

	# make egene dictionaries for each file 

	gtex_egenes = get_lead_snps(gtex_allpairs, 6, 7)
	admixed_local_egenes = get_lead_snps(allanc_allpairs, 8, 9)
	admixed_global_egenes = get_lead_snps(allanc_allpairs, 5, 6)
	admixed_allanc_egenes = get_lead_snps(allanc_allpairs, 14, 15)
	
	if allsubj_allpairs is not False:
		allsubj_local_egenes = get_lead_snps(allsubj_allpairs, 8, 9)
		allsubj_global_egenes = get_lead_snps(allsubj_allpairs, 5, 6) # should be the same as GTEx v8 results (gtex_egenes)
	else:
		allsubj_local_egenes = False
		allsubj_global_egenes = False 

	# save dictionaries 

	outdir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons'
	if tissue == 'Skin_Not_Sun_Exposed_Suprapubic':
		dicts = [ gtex_egenes, admixed_local_egenes, admixed_global_egenes, admixed_allanc_egenes, allsubj_local_egenes, allsubj_global_egenes ]
		names = [ '/gtex_v8/gtex_egenes', '/both_anc/admixed_local_egenes', '/both_anc/admixed_global_egenes', '/both_anc/admixed_allanc_egenes', '/all_ind/allsubj_local_egenes', '/all_ind/allsubj_global_egenes' ]
	else: 
		dicts = [ gtex_egenes, admixed_local_egenes, admixed_global_egenes, admixed_allanc_egenes ]
		names = [ '/gtex_v8/gtex_egenes', '/both_anc/admixed_local_egenes', '/both_anc/admixed_global_egenes', '/both_anc/admixed_allanc_egenes' ]	

	for i in range(len(dicts)):
		with open(outdir + names[i] + '.' + tissue + '.pkl', 'wb') as handle:
			pickle.dump(dicts[i], handle, protocol=pickle.HIGHEST_PROTOCOL)

	if tissue == 'Muscle_Skeletal':
		return [gtex_egenes, admixed_local_egenes, admixed_global_egenes, admixed_allanc_egenes]
	else:
		return [gtex_egenes, admixed_local_egenes, admixed_global_egenes, admixed_allanc_egenes, allsubj_local_egenes, allsubj_global_egenes]



def load_pickle(tissue):
	
	outdir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons'
	if tissue == 'Muscle_Skeletal':
		names = [ '/gtex_v8/gtex_egenes', '/both_anc/admixed_local_egenes', '/both_anc/admixed_global_egenes', '/both_anc/admixed_allanc_egenes' ]
	else:
		names = [ '/gtex_v8/gtex_egenes', '/both_anc/admixed_local_egenes', '/both_anc/admixed_global_egenes', '/both_anc/admixed_allanc_egenes', '/all_ind/allsubj_local_egenes', '/all_ind/allsubj_global_egenes' ]

	pickle_list = [0] * len(names)
	for i in range(len(names)):
		with open(outdir + names[i] + '.' + tissue + '.pkl', 'rb') as handle:
			b = pickle.load(handle)
		pickle_list[i] = b
	
	return pickle_list
	

if __name__ == '__main__':

	# merged = {}

	# tissue = sys.argv[1]

	# if os.path.isfile('/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/gtex_v8/gtex_egenes' + '.' + tissue + '.pkl'):
	# 	if tissue == 'Muscle_Skeletal':
	# 		gtex_egenes, admixed_local, admixed_global, admixed_allanc = load_pickle(tissue)
	# 	else:
	# 		gtex_egenes, admixed_local, admixed_global, admixed_allanc, allsubj_local, allsubj_global, newpeers_local, newpeers_global = load_pickle(tissue)	

	# else:
	# 	if tissue == 'Muscle_Skeletal':
	# 		gtex_egenes, admixed_local, admixed_global, admixed_allanc = parse_files(tissue)
	# 	else:
	# 		gtex_egenes, admixed_local, admixed_global, admixed_allanc, allsubj_local, allsubj_global, newpeers_local, newpeers_global = parse_files(tissue)

	# # rename slope and pvalue columns 
	# if tissue == "Skin_Not_Sun_Exposed_Suprapubic":
	# 	newcol = [ '_gtex', '_admixed_local', '_admixed_global', '_admixed_allanc', '_allsubj_local', '_allsubj_global' ]
	# 	dictlist = [gtex_egenes, admixed_local, admixed_global, admixed_allanc, allsubj_local, allsubj_global]
	# else:
	# 	newcol = [ '_gtex', '_admixed_local', '_admixed_global', '_admixed_allanc' ]
	# 	dictlist = [gtex_egenes, admixed_local, admixed_global, admixed_allanc]		
	
	# dflist = [0] * len(dictlist)
	# for i in range(len(dictlist)):
	# 	# convert dict to df 
	# 	df = dict_to_long_data(dictlist[i])
	# 	dflist[i] = df.rename(index=str, columns={'pval_nominal': 'pval_nominal'+newcol[i], 'slope': 'slope'+newcol[i], 'count': 'count'+newcol[i]})

	# if tissue == "Skin_Not_Sun_Exposed_Suprapubic":
	# 	gtex_egenes, admixed_local, admixed_global, admixed_allanc, allsubj_local, allsubj_global = dflist
	# else:
	# 	gtex_egenes, admixed_local, admixed_global, admixed_allanc = dflist

	# # merge data frames
	# m = gtex_egenes.merge(admixed_local, how='outer', on=['gene_id','variant_id'])
	# m = m.merge(admixed_global, how='outer', on=['gene_id','variant_id'])
	# m = m.merge(admixed_allanc, how='outer', on=['gene_id','variant_id'])
	# if tissue == 'Skin_Not_Sun_Exposed_Suprapubic':
	# 	m = m.merge(allsubj_local, how='outer', on=['gene_id','variant_id'])
	# 	m = m.merge(allsubj_global, how='outer', on=['gene_id','variant_id'])

	# m.to_csv('/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/chr19-merged-all-comparisons-tied-egenes-' +tissue + '.tsv.gz', sep='\t', header=True, index=False, compression='gzip')

	# print m.head()



	# compare only original version to newpeers version in Skin
	merged = {}
	tissue = 'Skin_Not_Sun_Exposed_Suprapubic'
	admixed_local, admixed_global, newpeers_local, newpeers_global = parse_newpeers(tissue)

	# rename slope and pvalue columns 
	newcol = [ '_admixed_local', '_admixed_global', '_newpeers_local', '_newpeers_global' ]
	dictlist = [admixed_local, admixed_global, newpeers_local, newpeers_global]
	
	dflist = [0] * len(dictlist)
	for i in range(len(dictlist)):
		# convert dict to df 
		df = dict_to_long_data(dictlist[i])
		dflist[i] = df.rename(index=str, columns={'pval_nominal': 'pval_nominal'+newcol[i], 'slope': 'slope'+newcol[i], 'count': 'count'+newcol[i]})

	admixed_local, admixed_global, newpeers_local, newpeers_global = dflist

	# merge data frames
	m = admixed_local.merge(admixed_global, how='outer', on=['gene_id','variant_id'])
	m = m.merge(newpeers_local, how='outer', on=['gene_id','variant_id'])
	m = m.merge(newpeers_global, how='outer', on=['gene_id','variant_id'])

	m.to_csv('/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr19_comparisons/chr19-merged-newpeers-tied-egenes-' +tissue + '.tsv.gz', sep='\t', header=True, index=False, compression='gzip')

	print m.head()

