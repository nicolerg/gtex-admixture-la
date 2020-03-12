# Filter allpairs file

import sys
import gzip

infile = sys.argv[1]
outfile = infile.replace('-allpairs','-allpairs-filt')

with gzip.open(infile, 'rb') as allpairs, gzip.open(outfile, 'wb') as out:
	out.write('gene_id\tvariant_id\ttss_distance\tma_samples\tma_count\tglobal_pval_nominal\tglobal_slope\tglobal_slope_se\tlava_pval_nominal\tlava_slope\tlava_slope_se\tafr_pval_nominal\tafr_slope\tafr_slope_se\n')
	next(allpairs)
	for line in allpairs:
		l = line.strip().split()
		try:
			global_pval = float(l[5])
			lava_pval = float(l[8])
			if lava_pval < 0.05 or global_pval < 0.05:
				out.write(line.strip()+'\n')
		except ValueError:
			continue

	