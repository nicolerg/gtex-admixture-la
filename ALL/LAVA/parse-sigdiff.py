import gzip
import os
import sys 

mypairs=sys.argv[1]
allpairs=sys.argv[2]
sigpairs=sys.argv[3]
outlocal=sys.argv[4]
outgtex=sys.argv[5]

# mypairs="/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/wholeblood-local-allpairs-chr21.tsv.gz"
# allpairs="/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/Whole_Blood.allpairs.chr21.txt.gz"
# sigpairs="/users/nicolerg/gtex-admix/fixed/wholeblood-chr21-sigdiff-pairs.tsv"
# outlocal="/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/wholeblood-chr21-local-sigdifflabels.tsv"
# outgtex="/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood/wholeblood-chr21-gtex-sigdifflabels.tsv"

genedict={}

with open(sigpairs,'r') as pairs:
    for line in pairs:
        line = line.strip().split('\t')
        if line[0] not in genedict:
        	genedict[line[0]] = []
        genedict[line[0]].append((line[0],line[1]))

with open(outgtex, 'w') as outfile:
	outfile.write("gene_id\tvariant_id\ttss_distance\tma_samples\tma_count\tmaf\tpval_nominal\tslope\tslope_se\tsigdiff\n")
	with gzip.open(allpairs,'r') as pairs:
		next(pairs)
		for l in pairs:
			line = l.strip().split('\t')
			if line[0] in genedict:
				if (line[0], line[1]) in genedict[line[0]]:
					outfile.write(l.strip()+"\t1\n")
				else:
					outfile.write(l.strip()+"\t0\n")
			else:
				outfile.write(l.strip()+"\t0\n")

with open(outlocal, 'w') as out:
	out.write("gene_id\tvariant_id\ttss_distance\tma_samples\tma_count\tmaf\tpval_nominal\tslope\tslope_se\tafr_pval_nominal\tafr_slope\tafr_slope_se\tasn_pval_nominal\tasn_slope\tasn_slope_se\tsigdiff\n")
	with gzip.open(mypairs,'r') as my:
		next(my)
		for l in my:
			line = l.strip().split('\t')
			if line[0] in genedict:
				if (line[0], line[1]) in genedict[line[0]]:
					out.write(l.strip()+"\t1\n")
				else:
					out.write(l.strip()+"\t0\n")
			else:
				out.write(l.strip()+"\t0\n")
