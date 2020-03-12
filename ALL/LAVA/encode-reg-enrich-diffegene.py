import gzip
import pandas as pd 
import sys
import os.path

"""
Annotate tied lead SNPs with ENCODE cis-Regulatory Elements 
Only include SNPs where tied SNPs are non-overlapping
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

# make a dictionary of all cRE annotations
with gzip.open(cre,'rb') as annot:
	for line in annot:
		l = line.strip().split('\t')
		chrom = l[0]
		start = int(l[1]) + 1
		stop = int(l[2]) + 1

		if chrom in reg:
			reg[chrom][start] = stop

outfile = outdir + '/' + tissue + '.diff.egene.annotated.tsv'
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

		annotation = 0
		if chrom in reg:
			for start in reg[chrom]:
				if start <= bp:
					if reg[chrom][start] >= bp:
						annotation = 1 				
						continue
				else:
					continue

		out.write(line.strip() + '\t' + str(annotation) + '\n')

