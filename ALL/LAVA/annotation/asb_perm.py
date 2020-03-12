import sys 
import gzip 
import pandas
import cPickle as pickle
import os.path

# permutation test for ASB 

annotation_dir = sys.argv[1]
cutoff = sys.argv[2]
tissue = sys.argv[3]

iterations=5000
summary = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/summary_'+cutoff+'.txt'
out = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/summary_'+tissue+'_'+cutoff+'_'+str(iterations)+'.pval'

df = pickle.load(open(annotation_dir+'/allele-specific-binding/gtex_mac10_snp_to_asb_map.pkl', 'rb'))
df.columns = ['rsid']


'''
# summary.txt
method  hit     N       NA
global  0       182     Adipose_Subcutaneous
LAVA    0       217     Adipose_Subcutaneous
global  1       2       Adipose_Subcutaneous
LAVA    1       1       Adipose_Subcutaneous
global  0       278     Nerve_Tibial
LAVA    0       264     Nerve_Tibial
global  1       1       Nerve_Tibial
LAVA    1       3       Nerve_Tibial
'''

n_actual_hits = {}
with open(summary, 'rb') as summ:
	next(summ)
	for line in summ:
		l = line.strip().split()
		if l[3] not in n_actual_hits:
			n_actual_hits[l[3]] = {}
		
		if l[0] == 'global' and l[1] == '0':
			n_actual_hits[l[3]]['global_tests'] = int(l[2])
 
		elif l[0] == 'global' and l[1] == '1':
			n_actual_hits[l[3]]['global'] = int(l[2])

		elif l[0] == 'LAVA' and l[1] == '0':
			n_actual_hits[l[3]]['local_tests'] = int(l[2])

		else:
			n_actual_hits[l[3]]['local'] = int(l[2])

'''
# need to replace this with reading in "summary"

n_actual_hits = {'Adipose_Subcutaneous' : {'global': 2,'local': 7, 'global_tests': 253, 'local_tests': 445},
	'Nerve_Tibial' : {'global': 3,'local': 6, 'global_tests': 456, 'local_tests': 405},
	'Artery_Tibial' : {'global': 3,'local': 7, 'global_tests': 303, 'local_tests': 806},
	'Muscle_Skeletal' : {'global': 1,'local': 7, 'global_tests': 269, 'local_tests': 347},
	'Lung' : {'global': 2,'local': 2, 'global_tests': 133, 'local_tests': 184},
	'Skin_Not_Sun_Exposed_Suprapubic' : {'global': 3,'local': 7, 'global_tests': 401, 'local_tests': 508}}
'''

with open(out, 'wb') as output: 

	# global

	if os.path.exists(annotation_dir + '/allele-specific-binding/'+tissue+'_global_perm_test_'+cutoff+'_'+str(iterations)+'.txt'):
		n_hits = [0] * iterations
		with open(annotation_dir + '/allele-specific-binding/'+tissue+'_global_perm_test_'+cutoff+'_'+str(iterations)+'.txt', 'rb') as hits:
			i = 0 
			for line in hits:
				l = line.strip()
				n_hits[i] = int(l)
				i += 1
	
	else:
		n_perm = n_actual_hits[tissue]['global_tests']
		n_hits = [0] * iterations
		with open(annotation_dir +'/allele-specific-binding/'+tissue+'_global_perm_test_'+cutoff+'_'+str(iterations)+'.txt','wb') as perm:
			for i in range(iterations):
				sub = df.sample(n=n_perm, axis=0)
				n_hit = len(sub.index[sub['rsid'] != 'NA'].tolist())
				n_hits[i] = n_hit
				print n_hit
				perm.write(str(n_hit)+'\n')
		
	n_hits_ordered = n_hits.sort()
	global_hits = n_actual_hits[tissue]['global']
	global_pval = 0
	for i in range(iterations):
		if n_hits[i] >= global_hits:
			global_pval = float(iterations - i)/iterations
			break

	# local

	if os.path.exists(annotation_dir + '/allele-specific-binding/'+tissue+'_local_perm_test_'+cutoff+'_'+str(iterations)+'.txt'):
		n_hits = [0] * iterations
		with open(annotation_dir + '/allele-specific-binding/'+tissue+'_local_perm_test_'+cutoff+'_'+str(iterations)+'.txt', 'rb') as hits:
			i = 0 
			for line in hits:
				l = line.strip('')
				n_hits[i] = int(l)
				i += 1

	else:
		n_perm = n_actual_hits[tissue]['local_tests']
		n_hits = [0] * iterations
		with open(annotation_dir + '/allele-specific-binding/'+tissue+'_local_perm_test_'+cutoff+'_'+str(iterations)+'.txt','wb') as perm:
			for i in range(iterations):
				sub = df.sample(n=n_perm, axis=0)
				n_hit = len(sub.index[sub['rsid'] != 'NA'].tolist())
				n_hits[i] = n_hit
				print n_hit
				perm.write(str(n_hit)+'\n')

	n_hits_ordered = n_hits.sort()
	local_hits = n_actual_hits[tissue]['local']
	local_pval = 0
	for i in range(iterations):
		if n_hits[i] >= local_hits:
			local_pval = float(iterations - i)/iterations
			break	

	output.write(tissue+'\n')
	output.write('global empirical p-value: '+str(global_pval)+'\n')
	output.write('local empirical p-value: '+str(local_pval)+'\n')
