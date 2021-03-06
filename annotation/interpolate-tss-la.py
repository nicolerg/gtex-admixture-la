# Interpolate local ancestry 
# Nicole Gay
import gzip
import sys

c = int(sys.argv[1])

tss = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/human_tss_hg19.bed'
mapfile = "/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/ancestry.chr"+str(c)+".bed"
outfile = "/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/chr"+str(c)+"_hg19_tss_localcov.tsv.gz"

# make a dictionary of dictionaries for all haplotypes
master = {}
with open(mapfile, 'r') as snpmap:
	next(snpmap)
	for line in snpmap:
		line = line.strip().split('\t')
		subjid = line[4].split('_')[0] # GTEx subject ID
		try:
			haplo = line[4].split('_')[1] # A or B
		except IndexError:
			print line
		start = line[1]
		stop = line[2]
		anc = line[3]
		if subjid in master:
			if haplo in master[subjid]:
				master[subjid][haplo][start] = [stop, anc]
			else:
				master[subjid][haplo] = {}
				master[subjid][haplo][start] = [stop, anc]
		else:
			master[subjid] = {}
			master[subjid][haplo] = {}
			master[subjid][haplo][start] = [stop, anc]

# now go through the subset file and assign local ancestry 
with open(tss, 'rb') as sub, gzip.open(outfile, 'wb') as out:
	# header
	out.write('SUBJID\tGENE_ID\tCHR\tPOS\tAFR\tASN\tEUR\tUNK\n')
	next(sub)
	# now find local ancestry for each SNP
	for line in sub:
		line = line.strip().split()
		chrom = line[0]
		if chrom == 'chr'+ str(c):
			pos = int(line[2])
			g = line[3]
			for key in master: # for each subject 
				# initialize local covariates
				AFR_1 = 0
				EUR_1 = 0
				ASN_1 = 0
				UNK_1 = 0
				AFR_2 = 0
				EUR_2 = 0
				ASN_2 = 0
				UNK_2 = 0
				# look at haplotype A/1
				adict = master[key]['A']
				anc = ''
				for start in adict:
					if pos >= int(start) and pos <= int(adict[start][0]):
						anc = adict[start][1]
						break 
				if anc == 'AFR':
					AFR_1 = 1
				elif anc == 'EUR':
					EUR_1 = 1
				elif anc == 'ASN':
					ASN_1 = 1
				else:
					UNK_1 = 1
				# look at haplotype B/2
				bdict = master[key]['B']
				anc = ''
				for start in bdict:
					if pos >= int(start) and pos <= int(bdict[start][0]):
						anc = bdict[start][1]
						break 
				if anc == 'AFR':
					AFR_2 = 1
				elif anc == 'EUR':
					EUR_2 = 1
				elif anc == 'ASN':
					ASN_2 = 1
				else:
					UNK_2 = 1
				AFR = AFR_1 + AFR_2
				EUR = EUR_1 + EUR_2
				ASN = ASN_1 + ASN_2
				UNK = UNK_1 + UNK_2
				out.write('\t'.join([key, g, chrom, str(pos), str(AFR), str(ASN), str(EUR), str(UNK)]) + '\n')
