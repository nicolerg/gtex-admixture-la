import gzip
import pandas as pd 
import sys
import os.path

"""
Annotate tied lead SNPs with ENCODE cis-Regulatory Elements 
Perform annotation separately for global and LAVA lead SNPs
Keep track of number of annotations per gene
"""

infile = sys.argv[1]
cre = sys.argv[2]
outdir = sys.argv[3]
tissue = sys.argv[4]
method = sys.argv[5]

## columns of infile:
# gene_id
# count
# variant_id
# pval_nominal
# slope
# gtex_cutoff

## columns of cRE annotation:
# chr, start, stop (0-indexed)

# initialize dictionary
reg = {}
for c in range(1,23):
	reg['chr'+str(c)] = {}

# initialize dictionary for annotation counts
gene_annot = {}

# make a dictionary of all cRE annotations
with gzip.open(cre,'rb') as annot:
	for line in annot:
		l = line.strip().split('\t')
		chrom = l[0]
		start = int(l[1]) + 1
		stop = int(l[2]) + 1

		if chrom in reg:
			reg[chrom][start] = stop

outfile = outdir + '/' + tissue + '.' + method + '.eSNP.annotated.tsv'
# iterate through input file; add boolean annotation column
with open(infile, 'rb') as ties, open(outfile, 'wb') as out:
	out.write('gene_id\tcount\tvariant_id\tpval_nominal\tslope\tgtex_cutoff\tannotation\n')
	next(ties)
	for line in ties:
		l = line.strip().split()
		variant_id = l[2]
		chrom = variant_id.split('_')[0]
		bp = int(variant_id.split('_')[1])

		gene_id = l[0]
		if gene_id in gene_annot:
			gene_annot[gene_id][1] += 1
		else:
			gene_annot[gene_id] = [0,1,l[3]]

		annotation = 0
		if chrom in reg:
			for start in reg[chrom]:
				if start <= bp:
					if reg[chrom][start] >= bp:
						annotation = 1 
						gene_annot[gene_id][0] += 1
						continue
				else:
					continue

		out.write(line.strip() + '\t' + str(annotation) + '\n')

# write out empirical probabilities of having an annotation in a top SNP for each gene
genes = [0] * len(gene_annot) 
probs = [0] * len(gene_annot) 
pval = [0] * len(gene_annot)
i = 0
for gene in gene_annot:
	genes[i] = gene
	total_count = gene_annot[gene][1]
	total_annot = gene_annot[gene][0]
	pval[i] = gene_annot[gene][2]
	if total_annot == 0:
		probs[i] = 0.0 
	else:
		probs[i] = float(total_annot)/float(total_count)
	i += 1

outfile = outdir + '/' + tissue + '.' + method + '.annot.prob.tsv'
df = pd.DataFrame({'gene_id':genes, 'prob_annot':probs, 'pval_nominal':pval})
df.to_csv(outfile, sep='\t', index=False)
