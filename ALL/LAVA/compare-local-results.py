# Author: Nicole Gay
# Winter 2018
# Input: files from slow_eqtl.R 
# Output: tables of nominal p-values for matching tests between LAVA and GTEx

import gzip
import os
import sys
import pandas as pd 

eqtldir = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/"
gtex_allpairs = eqtldir + "GTEx_Analysis_v8_eQTL_all_associations/Lung.allpairs.txt.gz"
gtex_eqtl = eqtldir + "GTEx_Analysis_v8_eQTL/Lung.v8.egenes.txt.gz"
gtex_allsig = eqtldir + "GTEx_Analysis_v8_eQTL/Lung.v8.signif_variant_gene_pairs.txt.gz"
indir = "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/"
my_allpairs = indir + "lung-LAVA-global-allpairs-merged.tsv.gz"
global_eqtl = indir + "global.eqtl.genes.tsv"
lava_eqtl = indir + "lava.eqtl.genes.tsv"
global_egenes = indir + "gtex.admix.global.egenes.txt.gz"

allsigout = indir + "gtex.allsig.global.allpairs.match.tsv.gz"

genedict={}

# compare my_allpairs to gtex_allsig

# make dictionary of gtex_allsig

with gzip.open(gtex_allsig, 'rb') as allsig:
    next(allsig)
    for line in allsig:
        l = line.strip().split()
        gene,variant,pval_nominal,effect = [ l[1], l[0], float(l[6]), float(l[7]) ]
        if gene not in genedict:
            genedict[gene] = {}
        genedict[gene][variant] = [ pval_nominal, effect ]

def check_val(val):
    if val != "ERROR" and val != "NA":
        return float(val)
    return None

# now iterate through my_allpairs and write out file of matches
with gzip.open(my_allpairs, 'rb') as allpairs, gzip.open(allsigout, 'wb') as out:
    out.write('\t'.join([ 'gene_id', 'variant_id', 'gtex_pval_nominal', 'gtex_slope', 'global_pval_nominal', 'global_slope', 'lava_pval_nominal', 'lava_slope' ]) + '\n')
    next(allpairs)
    for line in allpairs: 
        l = line.strip().split()
        gene,variant,g_pval,g_slope,lava_pval,lava_slope = [ l[0], l[1], l[5], l[6], l[8], l[9] ]
        lava_pval = check_val(lava_pval)
        lava_slope = check_val(lava_slope)
        g_pval = check_val(g_pval)
        g_slope = check_val(g_slope)
        if lava_pval is None or lava_slope is None or g_pval is None or g_slope is None:
            continue 
        if gene in genedict:
            if variant in genedict[gene]:
                gtex_pval = genedict[gene][variant][0]
                gtex_slope = genedict[gene][variant][1]
                out.write('\t'.join([gene,variant,str(gtex_pval),str(gtex_slope),str(g_pval),str(g_slope),str(lava_pval),str(lava_slope)]) + '\n')


# compare global eQTLs and LAVA eQTLs to GTEx allpairs 

# make dictionaries of LAVA and global eqtls

def make_dict(file):
    d = {}
    with open(file, 'rb') as eqtl:
        next(eqtl)
        for line in eqtl:
            l = line.strip().split()
            gene,variant,pval,slope = [ l[0], l[5], float(l[2]), float(l[7]) ]
            d[gene] = {}
            d[gene]['variant'] = variant
            d[gene]['pval_nominal'] = pval
            d[gene]['slope'] = slope
    return d

lavad = make_dict(lava_eqtl)
globald = make_dict(global_eqtl)

# iterate through GTEx allpairs

def check_dict(d,line):
    l = line.strip().split()
    gene,variant,pval,s = [ l[0], l[1], float(l[6]), float(l[7])]
    if gene in d:
        myvar = d[gene]['variant']
        if myvar == variant:
            d[gene]['gtex_pval_nominal'] = pval
            d[gene]['gtex_slope'] = s
    return d

with gzip.open(gtex_allpairs, 'rb') as allp:
    next(allp)
    for line in allp:
        lavad = check_dict(lavad,line)
        globald = check_dict(globald,line)

# write out file
gg = pd.DataFrame.from_dict(globald,orient='index')
gg.to_csv(indir+'global.eqtl.gtex.allpairs.compare.txt.gz',sep='\t',index=True,index_label='gene_id',compression='gzip')

lg = pd.DataFrame.from_dict(lavad,orient='index')
lg.to_csv(indir+'lava.eqtl.gtex.allpairs.compare.txt.gz',sep='\t',index=True,index_label='gene_id',compression='gzip')

