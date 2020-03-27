import sys

infile = 'gwas_char.txt'
config = 'config.json'

with open(infile, 'rb') as meta, open(config, 'wb') as out:
	next(meta)
	for line in meta:
		l = line.strip().split()
		# file path 
		if not l[1].startswith('WojcikG'):
			file = 'coloc_imputed_{}.txt.gz'.format(l[1])
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
