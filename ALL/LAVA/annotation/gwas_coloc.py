# get GWAS colocs from unique eGenes

import sys 
import gzip 
import pandas 
import os 
import pickle

# # curated from output of get_uniq_egenes.R; any uniq eGenes at 1e-5 or 1e-6
# base = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/master'
# gene_list = base + '/master.uniq.egene.txt'
# list_1e_5 = base + '/uniq_egenes_per_tissue_1e-05.tsv'
# list_1e_6 = base + '/uniq_egenes_per_tissue_1e-06.tsv'

# coloc_dir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/results_enloc_priors'

# tissues = ['Adipose_Subcutaneous', 'Nerve_Tibial', 'Artery_Tibial', 'Muscle_Skeletal', 'Lung', 'Skin_Not_Sun_Exposed_Suprapubic']

# master = {}
# for t in tissues:
# 	master[t] = {}

# gene_set = set()
# with open(gene_list, 'rb') as genes:
# 	for line in genes:
# 		l = line.strip().split()
# 		gene_set.add(l[0])

# # iterate through GWAS coloc files 

# for gwas in os.listdir(coloc_dir):
	
# 	tissue = gwas.strip().split('__')[2].split('.')[0]
# 	if tissue not in master:
# 		print 'Tissue ' + tissue + ' not found'
# 		continue

# 	phenotype = gwas.strip().split('__')[0]
# 	print phenotype

# 	infile = coloc_dir + '/' + gwas 

# 	with gzip.open(infile, 'rb') as coloc:
# 		next(coloc)
# 		for line in coloc:
# 			l = line.strip().split()
# 			gene_id = l[0]
# 			if gene_id in gene_set:
# 				if gene_id not in master[tissue]:
# 					master[tissue][gene_id] = {}
# 				try:
# 					master[tissue][gene_id][phenotype] = float(l[5]) # p4 
# 				except ValueError:
# 					master[tissue][gene_id][phenotype] = None # p4 

# # H0: neither trait has a genetic association in the region
# # H1: only trait 1 has a genetic association in the region
# # H2: only trait 2 has a genetic association in the region
# # H3: both traits are associated, but with different causal variants
# # H4: both traits are associated and share a single causal variant

# # gene_id p0      p1      p2      p3      p4
# # ENSG00000227232.5       0.009115167052395881    0.0003844977896782349   0.9504178524804027      0.0400
# # ENSG00000268903.1       0.1937050041906924      0.014112055556089822    0.7383982099410409      0.0537

# with open('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_egene_to_gwas_map.pkl', 'wb') as handle:
# 	pickle.dump(master, handle, protocol=pickle.HIGHEST_PROTOCOL)

master = pickle.load( open( '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_egene_to_gwas_map.pkl', 'rb' ) )

for key in master:
	outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_egene_to_gwas_map_'+key+'.tsv'
	d = master[key]
	df = pandas.DataFrame.from_dict(d, orient='index')
	df['tissue'] = key
	df.to_csv(outfile, sep='\t', index=True, index_label='gene_id', header=True)
