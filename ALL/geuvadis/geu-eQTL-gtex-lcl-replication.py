import gzip

indir = "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl"

highafrfile = indir+"/GEU-eQTLS-highAFR-effect.tsv"
higheurfile = indir+"/GEU-eQTLS-highEUR-effect.tsv"

outfile = indir+"/GEU-GTEX-LCL-diff-rvalue-plotdf.tsv"

# make dictionaries of the gene-SNP pairs for each group
highafr = {}
with open(highafrfile, 'rb') as afr:
	next(afr)
	for line in afr:
		l = line.strip().split()
		genebase = l[0].split('.')[0]
		eur_snp = l[12]
		afr_snp = l[24]
		highafr[genebase] = {}
		highafr[genebase]["EUR_ID"] = eur_snp
		highafr[genebase]["AFR_ID"] = afr_snp

higheur = {}
with open(higheurfile, 'rb') as eur:
	next(eur)
	for line in eur:
		l = line.strip().split()
		genebase = l[0].split('.')[0]
		eur_snp = l[12]
		afr_snp = l[24]
		higheur[genebase] = {}
		higheur[genebase]["EUR_ID"] = eur_snp
		higheur[genebase]["AFR_ID"] = afr_snp

# iterate through LCL allpairs

highafr_eursnp = 0
highafr_afrsnp = 0
higheur_eursnp = 0
higheur_afrsnp = 0

infile = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz"
with gzip.open(infile,'rb') as ap:
	next(ap)
	for line in ap:
		l = line.strip().split()
		if float(l[6]) < 0.05:
			genebase = l[0].split('.')[0]
			hg38_id = ':'.join(l[1].split('_')[0:2])
			if genebase in highafr:
				if highafr[genebase]["EUR_ID"] == hg38_id:
					highafr_eursnp += 1
					print highafr_eursnp
				if highafr[genebase]["AFR_ID"] == hg38_id:
					highafr_afrsnp += 1
					print highafr_afrsnp
			if genebase in higheur:
				if higheur[genebase]["EUR_ID"] == hg38_id:
					higheur_eursnp += 1
					print higheur_eursnp
				if higheur[genebase]["AFR_ID"] == hg38_id:
					higheur_afrsnp += 1
					print higheur_afrsnp

highafr_eursnp = float(highafr_eursnp)/float(len(highafr))
highafr_afrsnp = float(highafr_afrsnp)/float(len(highafr))
higheur_eursnp = float(higheur_eursnp)/float(len(higheur))
higheur_afrsnp = float(higheur_afrsnp)/float(len(higheur))

print 'Higher AFR effect size, AFR SNP-gene pair (fraction): ' + str(highafr_afrsnp)
print 'Higher AFR effect size, EUR SNP-gene pair (fraction): ' + str(highafr_eursnp)
print 'Higher EUR effect size, EUR SNP-gene pair (fraction): ' + str(higheur_eursnp)
print 'Higher EUR effect size, AFR SNP-gene pair (fraction): ' + str(higheur_afrsnp)

with open(outfile, 'wb') as out:
	out.write('POP\tGROUP\tFRAC\n')
	out.write('AFR\thighafr\t'+str(highafr_afrsnp)+'\n')
	out.write('EUR\thighafr\t'+str(highafr_eursnp)+'\n')
	out.write('AFR\thigheur\t'+str(higheur_afrsnp)+'\n')
	out.write('EUR\thigheur\t'+str(higheur_eursnp)+'\n')
