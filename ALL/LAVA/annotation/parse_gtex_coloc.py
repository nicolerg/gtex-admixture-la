# parse GTEx coloc results for significant colocalizations 

import os
import gzip 

indir='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_coloc/results_enloc_priors'
outfile='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_coloc/results_p4_gt_0.5.txt.gz'

with gzip.open(outfile, 'wb') as out:
	
	out.write('tissue\ttrait\tgene_id\tp4\n')

	for filename in os.listdir(indir):

		print filename

		full = os.path.join(indir, filename)
	
		# UKB_2395_2_Hair_or_balding_pattern_Pattern_2__PM__Liver.txt.gz
		trait = filename.split('__PM__')[0]
		tissue = filename.split('__PM__')[1].split('.')[0]

		with gzip.open(full, 'rb') as coloc:

			next(coloc)

			for line in coloc:
				l = line.strip().split()
				gene_id = l[0]
				p4 = l[5]
				
				if p4 == 'NA':
					continue

				if float(p4) > 0.5:
					out.write('\t'.join([tissue, trait, gene_id, p4]) + '\n')

