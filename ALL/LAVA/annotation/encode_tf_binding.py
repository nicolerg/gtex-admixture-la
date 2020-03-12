# annotate SNPs based on ENCODE TF-binding CHIP-seq
# https://www.encodeproject.org/matrix/?type=Experiment&assay_slims=DNA+binding&assay_title=ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&target.investigated_as=transcription+factor&biosample_ontology.classification=tissue&files.file_type=bed+narrowPeak

import gzip
import sys 

tissue = sys.argv[1]

labels = {'Adipose_Subcutaneous':'subcutaneous_adipose_tissue',
	'Muscle_Skeletal':'gastrocnemius_medialis',
	'Nerve_Tibial':'tibial_nerve',
	'Artery_Tibial':'tibial_artery',
	'Skin_Not_Sun_Exposed_Suprapubic':'suprapubic_skin',
	'Lung':'lung'}

peakdir = labels[tissue]

peaks = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/'+peakdir+'/'+peakdir+'.merged.bed.gz'
global_snps = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/'+tissue+'/gtex.admix.global.egenes.tied.txt.gz'
lava_snps = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/'+tissue+'/gtex.admix.lava.egenes.tied.txt.gz'

# gene_id count   variant_id      pval_nominal    slope   gtex_cutoff
# ENSG00000242268.2       4831    chr3_169273350_G_A_b38  0.00107171953022        1.12195926052   0.00014118

# initialize dictionary
reg = {}
for c in range(1,23):
	reg['chr'+str(c)] = {}

# make a dictionary of all peaks
with gzip.open(peaks,'rb') as tf:
	for line in tf:
		l = line.strip().split('\t')
		chrom = l[0]
		start = int(l[1]) + 1
		stop = int(l[2]) + 1

		if chrom in reg:
			reg[chrom][start] = stop

outfile = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/'+peakdir+'/'+peakdir+'.tf.binding.tsv'

# iterate through peak file; add boolean annotation column
with gzip.open(global_snps, 'rb') as gties, gzip.open(lava_snps, 'rb') as lties, open(outfile, 'wb') as out:
	out.write('gene_id\tcount\tvariant_id\tpval_nominal\tslope\tmethod\tannotation\n')

	next(gties)
	method='global'
	for line in gties:
		l = line.strip().split()
		variant_id = l[2]
		chrom = variant_id.split('_')[0]
		bp = int(variant_id.split('_')[1])

		annotation = 0
		if chrom in reg:
			for start in reg[chrom]:
				if start <= bp:
					if reg[chrom][start] >= bp:
						annotation = 1 				
						continue
				else:
					continue

		out.write(line.strip() + '\t' + method + '\t' + str(annotation) + '\n')

	next(lties)
	method='LAVA'
	for line in lties:
		l = line.strip().split()
		variant_id = l[2]
		chrom = variant_id.split('_')[0]
		bp = int(variant_id.split('_')[1])

		annotation = 0
		if chrom in reg:
			for start in reg[chrom]:
				if start <= bp:
					if reg[chrom][start] >= bp:
						annotation = 1 				
						continue
				else:
					continue

		out.write(line.strip() + '\t' + method + '\t' + str(annotation) + '\n')	

