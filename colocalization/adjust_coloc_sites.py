#!/bin/python3

import sys

# one file per gwas
infile = sys.argv[1]
outfile = infile.split('.formatted')[0] + '.adjusted.coloc-tests.txt'

# make sure the same lead SNP is used to test both eQTL files per tissue

# chr    snp_pos    gwas_file    eqtl_file    trait    gwas_pvalue    eqtl_pvalue    feature
# 17    39905964    /oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/gwas/imputed_TAGC_Asthma_EUR.formatted.txt.gz    /oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/allpairs_difflead/Adipose_Subcutaneous.filtered.global.allpairs.sorted.tsv.gz    /oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/gwas/imputed_TAGC_Asthma_EUR.formatted.txt.gz    1e-42    5.40886299272e-06    ENSG00000274512.5
# 17    39905964    /oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/gwas/imputed_TAGC_Asthma_EUR.formatted.txt.gz    /oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/allpairs_difflead/Adipose_Subcutaneous.filtered.global.allpairs.sorted.tsv.gz    /oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_COLOC/gwas/imputed_TAGC_Asthma_EUR.formatted.txt.gz    1e-42    1.06150278531e-05    ENSG00000274419.6

# first, read all of the hits into a dictionary

def get_tissue(file):
    for tissue in ["Lung", "Adipose_Subcutaneous", "Muscle_Skeletal", "Nerve_Tibial", "Whole_Blood", "Artery_Tibial", "Skin_Not_Sun_Exposed_Suprapubic"]:
        if str(tissue) in str(file):
            return tissue

def get_method(file):
    for method in ["global", "local"]:
        if str(method) in str(file):
            return method

seed_snps = {}
for tissue in ["Lung", "Adipose_Subcutaneous", "Muscle_Skeletal", "Nerve_Tibial", "Whole_Blood", "Artery_Tibial", "Skin_Not_Sun_Exposed_Suprapubic"]:
    seed_snps[tissue] = {}

with open(infile,'r') as loci:
    next(loci)
    for line in loci:
        l = line.strip().split()
        l = [str(i) for i in l]
        chrom, snp_pos, gwas_file, eqtl_file, trait, gwas_pvalue, eqtl_pvalue, feature = l
        tissue = get_tissue(eqtl_file)
        method = get_method(eqtl_file)
        if feature not in seed_snps[tissue]:
            seed_snps[tissue][feature] = {}
        if method not in seed_snps[tissue][feature]:
            seed_snps[tissue][feature][method] = {}
        seed_snps[tissue][feature][method]['chrom'] = chrom
        seed_snps[tissue][feature][method]['snp_pos'] = snp_pos
        seed_snps[tissue][feature][method]['gwas_file'] = str(gwas_file)
        seed_snps[tissue][feature][method]['eqtl_file'] = str(eqtl_file)
        seed_snps[tissue][feature][method]['gwas_pvalue'] = gwas_pvalue
        seed_snps[tissue][feature][method]['eqtl_pvalue'] = eqtl_pvalue

same_seed = 0
different_seed = 0
unique_loci = 0

# now iterate through the dictionary and output two lines for each locus
with open(outfile,'w') as out:
    out.write('chr\tsnp_pos\tgwas_file\teqtl_file\ttrait\tgwas_pvalue\teqtl_pvalue\tfeature\n')
    for tissue in seed_snps:
        for feature in seed_snps[tissue]:
            if len(seed_snps[tissue][feature]) == 2:
                # both methods say the locus should be tested
                # see if they're the same seed SNP
                global_seed = int(seed_snps[tissue][feature]['global']['snp_pos'])
                local_seed = int(seed_snps[tissue][feature]['local']['snp_pos'])
                if global_seed == local_seed:
                    same_seed += 1
                else:
                    different_seed += 1
                # different seed SNPs. test both
                for seed in set([global_seed, local_seed]): # either 1 or 2
                    # output line - global
                    chrom = seed_snps[tissue][feature]['global']['chrom']
                    gwas_file = seed_snps[tissue][feature]['global']['gwas_file']
                    gwas_pvalue = seed_snps[tissue][feature]['global']['gwas_pvalue']
                    eqtl_pvalue = seed_snps[tissue][feature]['global']['eqtl_pvalue']
                    eqtl_file = seed_snps[tissue][feature]['global']['eqtl_file']
                    out.write(('\t').join([chrom, str(seed), gwas_file, eqtl_file, gwas_file, gwas_pvalue, eqtl_file, feature]) + '\n')
                    # output line - local
                    chrom = seed_snps[tissue][feature]['local']['chrom']
                    gwas_file = seed_snps[tissue][feature]['local']['gwas_file']
                    gwas_pvalue = seed_snps[tissue][feature]['local']['gwas_pvalue']
                    eqtl_pvalue = seed_snps[tissue][feature]['local']['eqtl_pvalue']
                    eqtl_file = seed_snps[tissue][feature]['local']['eqtl_file']
                    out.write(('\t').join([chrom, str(seed), gwas_file, eqtl_file, gwas_file, gwas_pvalue, eqtl_file, feature]) + '\n')
            else:
                # only one method says the locus should be tested
                # test the same seed SNP for both methods
                unique_loci += 1
                if 'global' in seed_snps[tissue][feature]:
                    main_method = 'global'
                    global_eqtl = seed_snps[tissue][feature]['global']['eqtl_file']
                    local_eqtl = global_eqtl.replace('global', 'local')
                elif 'local' in seed_snps[tissue][feature]:
                    main_method = 'local'
                    local_eqtl = seed_snps[tissue][feature]['local']['eqtl_file']
                    global_eqtl = local_eqtl.replace('local', 'global')
                else:
                    print("I shouldn't get here.")

                chrom = seed_snps[tissue][feature][main_method]['chrom']
                snp_pos = seed_snps[tissue][feature][main_method]['snp_pos']
                gwas_file = seed_snps[tissue][feature][main_method]['gwas_file']
                gwas_pvalue = seed_snps[tissue][feature][main_method]['gwas_pvalue']
                eqtl_pvalue = seed_snps[tissue][feature][main_method]['eqtl_pvalue']
                # write line - global
                out.write(('\t').join([chrom, snp_pos, gwas_file, global_eqtl, gwas_file, gwas_pvalue, eqtl_file, feature]) + '\n')
                # write line - local
                out.write(('\t').join([chrom, snp_pos, gwas_file, local_eqtl, gwas_file, gwas_pvalue, eqtl_file, feature]) + '\n')

print('Same seed: ' + str(same_seed))
print('Different seed: ' + str(different_seed))
print('Unique loci: ' + str(unique_loci))
