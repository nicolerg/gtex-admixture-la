# get all GTEx v8 gene-level colocalization p-values

# DEPRECATED

import sys 
import gzip 
import pandas 
import os 
import pickle

# keep all genes 

coloc_dir = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/results_enloc_priors'

tissues = ['Adipose_Subcutaneous', 'Nerve_Tibial', 'Artery_Tibial', 'Muscle_Skeletal', 'Lung', 'Skin_Not_Sun_Exposed_Suprapubic']

outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master-coloc-v2.tsv'

with open(outfile, 'wb') as out:
	out.write('\t'.join(['gene_id','p3','p4','tissue','trait']) + '\n')

	for gwas in os.listdir(coloc_dir):
		
		tissue = gwas.strip().split('__')[2].split('.')[0]
		if tissue not in tissues:
			continue

		phenotype = gwas.strip().split('__')[0]
		print phenotype

		infile = coloc_dir + '/' + gwas 

		with gzip.open(infile, 'rb') as coloc:
			next(coloc)
			for line in coloc:
				l = line.strip().split()
				gene_id = l[0]
				try:
					p4 = float(l[5]) # p4 
				except ValueError:
					p4 = None # p4 
				try:
					p3 = float(l[4])
				except ValueError:
					p3 = None
				out.write('\t'.join([gene_id, str(p3), str(p4), tissue, phenotype]) + '\n')

# H0: neither trait has a genetic association in the region
# H1: only trait 1 has a genetic association in the region
# H2: only trait 2 has a genetic association in the region
# H3: both traits are associated, but with different causal variants
# H4: both traits are associated and share a single causal variant

# gene_id p0      p1      p2      p3      p4
# ENSG00000227232.5       0.009115167052395881    0.0003844977896782349   0.9504178524804027      0.0400
# ENSG00000268903.1       0.1937050041906924      0.014112055556089822    0.7383982099410409      0.0537

