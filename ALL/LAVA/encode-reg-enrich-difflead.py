import gzip
import pandas as pd 
import sys
import os.path

"""
Annotate tied lead SNPs with ENCODE cis-Regulatory Elements 
Only include SNPs where eGene is the same for both methods but tied SNPs are non-overlapping
Keep track of number of annotations per gene per method 

${tissue}.diff.SNP.annot.tsv --> gene-level summary of both methods
${tissue}.diff.eSNP.annotated.tsv --> one line per SNP per method per gene
"""

infile = sys.argv[1]
cre = sys.argv[2]
outdir = sys.argv[3]
tissue = sys.argv[4]

## columns of infile:
# gene_id
# count
# variant_id
# pval_nominal
# slope
# gtex_cutoff
# method

## columns of cRE annotation:
# chr, start, stop (0-indexed)

# initialize dictionary
reg = {}
for c in range(1,23):
	reg['chr'+str(c)] = {}

# initialize dictionaries for annotation counts
gene_annot = {} 

lava_annot = {}
global_annot = {}

# make a dictionary of all cRE annotations
with gzip.open(cre,'rb') as annot:
	for line in annot:
		l = line.strip().split('\t')
		chrom = l[0]
		start = int(l[1]) + 1
		stop = int(l[2]) + 1

		if chrom in reg:
			reg[chrom][start] = stop

outfile = outdir + '/' + tissue + '.diff.eSNP.annotated.tsv'
# iterate through input file; add boolean annotation column
with open(infile, 'rb') as ties, open(outfile, 'wb') as out:
	out.write('gene_id\tcount\tvariant_id\tpval_nominal\tslope\tgtex_cutoff\tmethod\tannotation\n')
	next(ties)
	for line in ties:
		l = line.strip().split()
		variant_id = l[2]
		chrom = variant_id.split('_')[0]
		bp = int(variant_id.split('_')[1])

		gene_id = l[0]
		
		if gene_id in gene_annot:
			gene_annot[gene_id][l[6]][2] = l[3]
			gene_annot[gene_id][l[6]][1] += 1 
		else:
			gene_annot[gene_id] = {}
			gene_annot[gene_id]['LAVA'] = [0,0,0]
			gene_annot[gene_id]['global'] = [0,0,0]

			gene_annot[gene_id][l[6]][2] = l[3] 
			gene_annot[gene_id][l[6]][1] += 1 

		annotation = 0
		if chrom in reg:
			for start in reg[chrom]:
				if start <= bp:
					if reg[chrom][start] >= bp:
						annotation = 1 
						gene_annot[gene_id][l[6]][0] += 1 					
						continue
				else:
					continue

		out.write(line.strip() + '\t' + str(annotation) + '\n')

# write out empirical probabilities of having an annotation in a top SNP for each gene
genes = [0] * len(gene_annot) 

lava_pval = [0] * len(gene_annot) 
global_pval = [0] * len(gene_annot) 

global_count = [0] * len(gene_annot) 
lava_count = [0] * len(gene_annot) 

global_annot = [0] * len(gene_annot) 
lava_annot = [0] * len(gene_annot) 

i = 0
for gene in gene_annot:
	genes[i] = gene 
	lava_pval[i] = gene_annot[gene]['LAVA'][2]
	global_pval[i] = gene_annot[gene]['global'][2]
	global_count[i] = gene_annot[gene]['global'][1]
	lava_count[i] = gene_annot[gene]['LAVA'][1]
	global_annot[i] = gene_annot[gene]['global'][0]
	lava_annot[i] = gene_annot[gene]['LAVA'][0]
	i += 1

outfile = outdir + '/' + tissue + '.diff.SNP.annot.tsv'
df = pd.DataFrame({'gene_id':genes, \
	'lava_pval':lava_pval, \
	'lava_count':lava_count, \
	'lava_annot':lava_annot, \
	'global_pval':global_pval, \
	'global_count':global_count, \
	'global_annot':global_annot})
df.to_csv(outfile, sep='\t', index=False)
