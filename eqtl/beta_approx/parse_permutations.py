infile = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/beta_permute.tsv'
outdir = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/chunks'
compare = '/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/genes_to_include_in_beta.tsv'

# make dictionary of number of tests/expected genes 
genes_to_test = {}
genes_to_test['Lung'] = {}
genes_to_test['Muscle_Skeletal'] = {}
genes_to_test['Adipose_Subcutaneous'] = {} # gene:num_tests

with open(compare, 'r') as comp:
	next(comp)
	for line in comp:
		l = line.strip().split()
		tissue = l[1]
		gene = l[0]
		n_snps = int(l[2])
		genes_to_test[tissue][gene] = {}
		genes_to_test[tissue][gene]['n_snps'] = n_snps
		genes_to_test[tissue][gene]['variants'] = set()

with open(infile, "r") as all_tests:  
    # Read next line
    header = all_tests.readline()
    # initialize iteration
    line = all_tests.readline()
    l = line.strip().split()
    current_tissue = l[0]
    current_gene = l[1]
    line_list = []
    snps_to_test = genes_to_test[current_tissue][current_gene]['n_snps']
    while line:
        l = line.strip().split()
        # tissue gene_id permutation variant_id global_pval_nominal local_pval_nominal
        gene = l[1]
        tissue = l[0]
        # do I want to test this gene?
        if gene in genes_to_test[tissue]:
        	# is this the current gene/tissue or a new one?
        	if gene == current_gene and tissue == current_tissue:
        		genes_to_test[tissue][gene]['variants'].add(l[3])
        		# don't append a line if we've already tested as many SNPs as we need to 
        		if len(genes_to_test[tissue][gene]['variants']) < snps_to_test:
        			line_list.append(line)
        		else:
        			print('{} {} {}'.format(current_tissue, current_gene, len(genes_to_test[tissue][gene]['variants'])))
        	else:
        		if(len(line_list) > 0):
        			#print('{} {} {}'.format(current_tissue, current_gene, len(line_list)))
        			with open('{}/{}_{}.txt'.format(outdir, current_tissue, current_gene), 'a') as f:
        				f.writelines(line_list)
        		# elif(len(line_list) > 0 and len(line_list) < 100):
        		# 	print('Why does {} {} have less than 100 lines?'.format(current_tissue, current_gene))
        		# 	print(line_list)
        		
        		# reinitialize 
        		current_tissue = tissue 
        		current_gene = gene 
        		snps_to_test = genes_to_test[tissue][gene]['n_snps']
        		line_list = []

        		line_list.append(line)
        		genes_to_test[tissue][gene]['variants'].add(l[3])
        		
        line = all_tests.readline()
