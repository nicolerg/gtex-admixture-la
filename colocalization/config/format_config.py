import sys

infile = 'gwas_char.txt'
all_gwas = 'gwas_experiments.json'
cc_only = 'gtex_cc.json'

with open(infile, 'rb') as meta, open(all_gwas, 'wb') as out, open(cc_only, 'wb') as cc:
	next(meta)
	for line in meta:
		l = line.strip().split()
		# file path 
		if not l[1].startswith('WojcikG'):
			file = 'imputed_{}.formatted.txt.gz'.format(l[1])
		else:
			file = l[1]

		out.write('"{}/{}":\n'.format(l[0], file))
		out.write('{\n')

		out.write('    "ref": "{}",\n'.format(l[5]))
		out.write('    "gwas_format": "{}",\n'.format(l[6]))
		
		gwas_type = l[3]
		N = l[2]
		out.write('    "N": {},\n'.format(N))
		if gwas_type == 'cc':
			s = l[4]
			out.write('    "s": {},\n'.format(s))
		out.write('    "type": "{}"\n'.format(gwas_type))
		out.write('},\n')

		## cc only:
		if not file.startswith('WojcikG'):

			if l[3] == 'cc':

				cc.write('"{}/{}":\n'.format(l[0], file))
				cc.write('{\n')
				cc.write('    "ref": "{}",\n'.format(l[5]))
				cc.write('    "gwas_format": "{}",\n'.format(l[6]))
				cc.write('    "N": {},\n'.format(l[2]))
				cc.write('    "s": {},\n'.format(l[4]))
				cc.write('    "type": "cc"\n')
				cc.write('},\n')

