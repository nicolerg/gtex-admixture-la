"""
Run this script on the merged allpairs file output by combined local & global eQTL
calling to summarize the top SNPs for each gene. All tied top eSNPs are output.

Currently, top SNPs for all tested genes are reported. They are not filtered by
significance. 
"""

import gzip
import pandas as pd 
import sys

infile = sys.argv[1]
tissue = sys.argv[2]
outdir = sys.argv[3]
# outdir = '/'.join( [ infile.strip().split('/')[:-1] ] ) + '/'

# gene_id
# variant_id
# # tss_distance
# ma_samples
# ma_count
# global_pval_nominal
# global_slope
# global_slope_se
# lava_pval_nominal
# lava_slope
# lava_slope_se
# afr_pval_nominal
# afr_slope
# afr_slope_se

# # make a dictionary of GTEx eGene cutoffs
# gtex_cutoffs = {}
# egenes = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/'+tissue+'.v8.egenes.txt.gz'
# with gzip.open(egenes,'rb') as e:
# 	next(e)
# 	for line in e:
# 		l = line.strip().split()
# 		gtex_cutoffs[l[0]] = float(l[29])

# make a dictionary of top SNPs 
global_ties = {}
local_ties = {}

def add_entry(d, line, pval_index, slope_index):
	l = line.strip().split()
	pval = l[pval_index]
	
	if l[0] in d:
		
		d[l[0]]['count'] += 1
		
		# break if pvalue is greater than current min 
		if float(pval) > d[l[0]]['pval_nominal']:
			return d

		# update dict if pvalue is less than current min
		if float(pval) < d[l[0]]['pval_nominal']:
			d[l[0]]['variant_id'] = [ l[1] ]
			d[l[0]]['pval_nominal'] = float(pval)
			d[l[0]]['slope'] = [ float(l[slope_index]) ]
			return d

		# append to existing entry if pvalues are the same
		if float(pval) == d[l[0]]['pval_nominal']:
			d[l[0]]['variant_id'].append(l[1])
			d[l[0]]['slope'].append(float(l[slope_index]))
			return d

	else: # new entry
		
		d[l[0]] = {}
		d[l[0]]['count'] = 1
		d[l[0]]['variant_id'] = [ l[1] ]
		d[l[0]]['pval_nominal'] = float(pval)
		d[l[0]]['slope'] = [ float(l[slope_index]) ]
		# d[l[0]]['tss_distance'] = int(float(l[2]))
		return d

fail = [ 'NaN', 'ERROR', 'NA']
with gzip.open(infile, 'rb') as allp:
	next(allp)
	for line in allp:
		l = line.strip().split()
		if l[5] in fail or l[8] in fail:
			continue
		global_ties = add_entry(global_ties, line, 5, 6)
		local_ties = add_entry(local_ties, line, 8, 9)

# # add the GTEx threshold
# for gene in global_ties:
# 	if gene in gtex_cutoffs:
# 		global_ties[gene]['gtex_cutoff'] = gtex_cutoffs[gene]
# 	else:
# 		global_ties[gene]['gtex_cutoff'] = 'NA'

# for gene in local_ties:
# 	if gene in gtex_cutoffs:
# 		local_ties[gene]['gtex_cutoff'] = gtex_cutoffs[gene]
# 	else:
# 		local_ties[gene]['gtex_cutoff'] = 'NA'

# write out results

def write_dict(d, outfile):
	tied_count = {}
	with gzip.open(outfile, 'wb') as out:
		out.write('gene_id\tcount\tvariant_id\tpval_nominal\tslope\n')
		for gene in d:
			gene_id = gene 
			pval = d[gene]['pval_nominal']
			c = d[gene]['count']
			tied_count[gene_id] = len(d[gene]['variant_id'])
			for i in range(len(d[gene]['variant_id'])):
				var = d[gene]['variant_id'][i]
				slope = d[gene]['slope'][i]
				out.write('\t'.join([ gene_id, str(c), var, str(pval), str(slope) ]) + '\n')
	return tied_count

global_counts = write_dict(global_ties, outdir+'/gtex.admix.global.egenes.tied.txt.gz')
local_counts = write_dict(local_ties, outdir+'/gtex.admix.local.egenes.tied.txt.gz')

# # make a data frame of counts
# gc = pd.DataFrame.from_dict(global_counts,orient='index',columns=['global_n_tied'])
# lc = pd.DataFrame.from_dict(local_counts,orient='index',columns=['lava_n_tied'])
# merged = gc.merge(lc, how='outer', left_index=True, right_index=True, copy=False)

# merged.to_csv(outdir+'/gtex.admix.tied.counts.txt',sep='\t',index=True,index_label='gene_id')

# gg = pd.DataFrame.from_dict(global_ties,orient='index')
# gg.to_csv(outdir+'gtex.admix.global.egenes.tied.txt.gz',sep='\t',index=True,index_label='gene_id',compression='gzip')

# lg = pd.DataFrame.from_dict(local_ties,orient='index')
# lg.to_csv(outdir+'gtex.admix.lava.egenes.tied.txt.gz',sep='\t',index=True,index_label='gene_id',compression='gzip')

# # now output eQTLs only 
# gg = gg[gg['pval_nominal'] <= gg['gtex_cutoff']]
# gg.to_csv(outdir+'gtex.admix.global.eqtl.txt.gz',sep='\t',index=True,index_label='gene_id',compression='gzip')

# lg = lg[lg['pval_nominal'] <= lg['gtex_cutoff']]
# lg.to_csv(outdir+'gtex.admix.lava.eqtl.txt.gz',sep='\t',index=True,index_label='gene_id',compression='gzip')

