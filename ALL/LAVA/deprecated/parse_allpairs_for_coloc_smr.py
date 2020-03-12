import gzip 
import sys 

tissue = sys.argv[1]
gene_list = sys.argv[2]
base = sys.argv[3]

allpairs = base + '/'+tissue+'/'+tissue+'-LAVA-global-allpairs-merged.tsv.gz'

af = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/snp_to_effect_af.tsv.gz'

local_out = base + '/annotation/coloc/allpairs_difflead/'+tissue+'.filtered.local.allpairs.tsv.gz'
global_out = base + '/annotation/coloc/allpairs_difflead/'+tissue+'.filtered.global.allpairs.tsv.gz'

genelist = set()
af_dict = {}

with gzip.open(af, 'rb') as snp_af:
	next(snp_af)
	for line in snp_af:
		l = line.strip().split()
		variant_id = l[0]
		af = float(l[3])
		af_dict[variant_id] = af

with open(gene_list, 'rb') as gl:
	for line in gl:
		g = line.strip()
		if g not in genelist:
			genelist.add(g)

with gzip.open(local_out, 'wb') as localp, gzip.open(global_out, 'wb') as globalp, gzip.open(allpairs, 'rb') as allp:

	localp.write('chr\tsnp_pos\talt\tref\tbeta\tse\tpvalue\tgene\teffect_af\n')
	globalp.write('chr\tsnp_pos\talt\tref\tbeta\tse\tpvalue\tgene\teffect_af\n')

	next(allp)

	for line in allp:

		l = line.strip().split()

		gene = l[0]

		if gene in genelist:

			chrom = l[1].strip().split('_')[0].replace('chr','')
			snp_pos = int(l[1].strip().split('_')[1])
			ref = l[1].strip().split('_')[2]
			alt = l[1].strip().split('_')[3]
			build = l[1].strip().split('_')[4]

			effect_af = af_dict[l[1]]
			
			try:
				pval_global = float(l[5])
				beta_global = float(l[6])
				se_global = float(l[7])
				
				pval_local = float(l[8])
				beta_local = float(l[9])
				se_local = float(l[10])

			except ValueError:
				continue 

			localp.write('\t'.join([chrom, str(snp_pos), alt, ref, str(beta_local), str(se_local), str(pval_local), gene, str(effect_af)]) + '\n')
			globalp.write('\t'.join([chrom, str(snp_pos), alt, ref, str(beta_global), str(se_global), str(pval_global), gene, str(effect_af)]) + '\n')
