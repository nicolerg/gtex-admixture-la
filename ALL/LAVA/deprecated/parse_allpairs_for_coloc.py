import gzip 
import sys 

tissue = sys.argv[1]

allpairs = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/'+tissue+'/'+tissue+'-LAVA-global-merged-allpairs.tsv.gz'
gene_list = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/'+tissue+'.1e-5.same.egene.diff.snp.txt'
local_out = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/'+tissue+'.filtered.local.allpairs.tsv.gz'
global_out = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/'+tissue+'.filtered.global.allpairs.tsv.gz'

genelist = set()

with open(gene_list, 'rb') as gl:
	for line in gl:
		g = line.strip()
		if g not in genelist:
			genelist.add(g)

with gzip.open(local_out, 'wb') as localp, gzip.open(global_out, 'wb') as globalp, gzip.open(allpairs, 'rb') as allp:

	localp.write('gene\tchr\tsnp_pos\tref\talt\tbuild\tbeta\tse\n')
	globalp.write('gene\tchr\tsnp_pos\tref\talt\tbuild\tbeta\tse\n')

	for line in allp:

		l = line.strip().split()

		gene = l[0]

		if gene in genelist:

			chrom = l[1].strip().split('_')[0].replace('chr','')
			snp_pos = int(l[1].strip().split('_')[1])
			ref = l[1].strip().split('_')[2]
			alt = l[1].strip().split('_')[3]
			build = l[1].strip().split('_')[4]
			
			try:
				beta_global = float(l[6])
			except ValueError:
				continue 

			try: 
				se_global = float(l[7])
			except ValueError:
				continue

			try:
				beta_local = float(l[9])
			except ValueError:
				continue

			try:
				se_local = float(l[10])
			except ValueError:
				continue 

			localp.write('\t'.join([gene, chrom, str(snp_pos), ref, alt, build, str(beta_local), str(se_local)]) + '\n')
			globalp.write('\t'.join([gene, chrom, str(snp_pos), ref, alt, build, str(beta_global), str(se_global)]) + '\n')
