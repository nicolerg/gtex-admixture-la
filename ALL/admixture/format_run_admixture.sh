# run ADMIXTURE on VCFs 
# following http://software.genetics.ucla.edu/admixture/admixture-manual.pdf

set -e

indir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814
outdir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/admixture_20190819

# ${indir} contains chromosome-level VCFs generated as follows:
# 1. GTEx and 1000 Genomes filtered with MAF > 0.05
# 2. find intersection of each chromosome
# 3. merge intersection
# 4. filter subject (117 admixed individuals and same reference panels used in RFMix) and apply MAF filter of 0.10
# 5. remove SNPs with missing values

# now do the following for ADMIXTURE:
# 1. LD-prune each chromosome. 
# 2. concatenate chromosomes
# 3. generate .pop file 
# 4. run ADMIXTURE 

# index 
for chr in {1..22}; do 
	tabix -p vcf ${indir}/merged_hg38_chr${chr}_filt_MAF10_nomissing.vcf.gz &
done
wait
# convert to binary .BED
for chr in {1..22}; do 
	plink --vcf ${indir}/merged_hg38_chr${chr}_filt_MAF10_nomissing.vcf.gz --make-bed --vcf-half-call m --out ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing & 
done
wait
# LD-prune each chromosome 
for chr in {1..22}; do 
	plink --bfile ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing --indep-pairwise 50 10 0.1 --out ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing &
done
wait
for chr in {1..22}; do 
	plink --bfile ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing --extract ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing.prune.in --make-bed --out ${outdir}/merged_hg38_chr${chr}_filt_MAF10_nomissing_pruned &
done
wait

cd ${outdir}

# concatenate all chromosomes 
ls | grep "pruned.bed" > bed
ls | grep "pruned.bim" > bim
ls | grep "pruned.fam" > fam
paste bed bim fam > files_to_merge.txt
sed -i -e '1d' files_to_merge.txt

plink --bfile merged_hg38_chr${chr}_filt_MAF10_nomissing_pruned --merge-list files_to_merge.txt --make-bed --out merged_hg38_filt_MAF10_nomissing_pruned

# remove intermediate files 
rm merged_hg38_chr*_filt_MAF10_nomissing.* merged_hg38_chr*_filt_MAF10_nomissing_pruned.* bed bim fam 

# 159817 variants 

# generate .pop file 
# I want to keep all individuals in this list: 
gtex_admixed=/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt
# as well as certain 1000 Genomes populations. I/we need to do multiple iterations with different reference panels, but include the following to start:
# Asian (N = 207)
#   CHB (Han Chinese in Beijing, China)
#   JPT (Japanese in Tokyo, Japan)
# African (N = 405)
#   YRI (Yoruba in Ibadan, Nigeria)
#   GWD (Gambian in Western Divisions in the Gambia)
#   MSL (Mende in Sierra Leone)
#   ESN (Esan in Nigeria)
# European (N = 99)
#   CEU
tg=/users/nicolerg/gtex-admix/metadata/1kg-subject-popcodes.tsv

Rscript ~/gtex-admix/scripts/admixture/order_pop_file.R ${outdir}/merged_hg38_filt_MAF10_nomissing_pruned.fam ${tg} 
cut -f3 merged_hg38_filt_MAF10_nomissing_pruned.pop.txt > merged_hg38_filt_MAF10_nomissing_pruned.pop

# add software to path 
export PATH=$PATH:/mnt/lab_data/montgomery/nicolerg/TOOLS/admixture_linux-1.3.0
# run admixture 
K=3 # number of reference populations included in the .pop file 
admixture merged_hg38_filt_MAF10_nomissing_pruned.bed ${K} --supervised

paste *pop.txt *.3.Q | sed "s/ /    /g" > merged_hg38_filt_MAF10_nomissing_pruned.3.Q.subj
cp merged_hg38_filt_MAF10_nomissing_pruned.3.Q.subj ~/gtex-admix/admixture/merged_hg38_filt_MAF10_nomissing_pruned_original-reference-panels.3.Q.subj

# plot interactively with plot_admixture.R

exit






















outdir=/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged # set to something in /srv/scratch or /mnt/lab_data/montgomery/<user>

# need to try a new workflow: first, merge ALL 1000 Genomes and GTEx SNPs
# ONLY THEN, LD-prune (I'm worried that thinning before mergind is leaving me with a poorly selected SNP set)
# do the merging and LD-pruning on a per-chromosome basis; run in parallel for each chromosome 

# 1000 Genomes VCFs, already split by chromosome:
tg=/mnt/lab_data/montgomery/shared/1KG/hg38

# GTEx VCF, needs to be split by chromosome:
gtex_vcf=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
# split GTEx VCF by chromosome (only need austosomes)
for chr in {1..22}; do 
	bcftools view -r chr${chr} -Oz -o ${outdir}/gtex_phased_chr${chr}.vcf.gz ${gtex_vcf} &
done
wait

# MAF-filter 
for chr in {1..22}; do 
	bcftools view --min-af 0.05 -Oz -o ${outdir}/gtex_phased_chr${chr}_MAF05.vcf.gz ${outdir}/gtex_phased_chr${chr}.vcf.gz &
	
done
wait
for chr in {1..22}; do 
	bcftools view --min-af 0.05 -Oz -o ${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz ${tg}/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz &
done
wait

# index VCFs 
for file in $(ls ${outdir} | grep -E "MAF05.vcf.gz$"); do 
	tabix -p vcf ${outdir}/${file} &
done

# merge each chromosome
for chr in {1..22}; do 
	tg_vcf=${outdir}/tg_phased_chr${chr}_MAF05.vcf.gz
	gtex_vcf=${outdir}/gtex_phased_chr${chr}_MAF05.vcf.gz
	bcftools merge -O z -o ${outdir}/merged_hg38_chr${chr}.vcf.gz ${tg_vcf} ${gtex_vcf} &
done
wait

# MAF-filter again
for chr in {1..22}; do 
	bcftools view --min-af 0.1 -Oz -o ${outdir}/merged_hg38_chr${chr}_MAF10.vcf.gz ${outdir}/merged_hg38_chr${chr}.vcf.gz &
done
wait

# # remove multiallelic variants
# for chr in {1..22}; do 
# 	bcftools view --max-alleles 2 -Oz -o ${outdir}/merged_hg38_chr${chr}_biallelic.vcf.gz ${outdir}/merged_hg38_chr${chr}.vcf.gz &
# done

# index
for chr in {1..22}; do 
	tabix -p vcf ${outdir}/merged_hg38_chr${chr}_MAF10.vcf.gz &
done
wait
# convert to binary .BED
for chr in {1..22}; do 
	plink --vcf ${outdir}/merged_hg38_chr${chr}_MAF10.vcf.gz --make-bed --vcf-half-call m --out ${outdir}/merged_hg38_chr${chr}_MAF10 & 
done

exit






# index
for file in $(ls ${outdir} | grep -E "MAF_0.1"); do tabix -p vcf ${outdir}/${file} & done
wait

# convert to plink 
for chr in {1..22}; do 
	plink --vcf ${outdir}/merged_hg38_chr${chr}_biallelic.vcf.gz --make-bed --vcf-half-call m --out ${outdir}/merged_hg38_chr${chr}_biallelic &
done
wait

# LD-prune each chromosome 
for chr in {1..22}; do 
	plink --bfile ${outdir}/merged_hg38_chr${chr}_biallelic --indep-pairwise 50 10 0.1 --out ${outdir}/merged_hg38_chr${chr}_biallelic &
done
wait
for chr in {1..22}; do 
	plink --bfile ${outdir}/merged_hg38_chr${chr}_biallelic --extract ${outdir}/merged_hg38_chr${chr}_biallelic.prune.in --make-bed --out ${outdir}/merged_hg38_chr${chr}_biallelic_pruned &
done
wait

# concatenate all chromosomes 
ls | grep "pruned.bed" > bed
ls | grep "pruned.bim" > bim
ls | grep "pruned.fam" > fam
paste bed bim fam > files_to_merge.txt
sed -i -e '1d' files_to_merge.txt

plink --bfile merged_hg38_chr11_biallelic_pruned --merge-list files_to_merge.txt --make-bed --out merged_hg38_biallelic_pruned

# remove missnps
for chr in {1..22}; do 
	plink --bfile merged_hg38_chr${chr}_biallelic_pruned_fixed --exclude merged_hg38_biallelic_pruned_fixed_pruned-merge.missnp --make-bed --out merged_hg38_chr${chr}_biallelic_pruned_fixed2 &
done

# concatenate all chromosomes 
ls | grep "pruned_fixed2.bed" > bed
ls | grep "pruned_fixed2.bim" > bim
ls | grep "pruned_fixed2.fam" > fam
paste bed bim fam > files_to_merge.txt
sed -i -e '1d' files_to_merge.txt

plink --bfile merged_hg38_chr10_biallelic_pruned_fixed2 --merge-list files_to_merge.txt --make-bed --out merged_hg38_biallelic_pruned_fixed_pruned

# LD-prune again 
plink --bfile merged_hg38_biallelic_pruned_fixed_pruned --indep-pairwise 70 10 0.1 --out merged_hg38_biallelic_pruned_fixed_pruned 
plink --bfile merged_hg38_biallelic_pruned_fixed_pruned --extract merged_hg38_biallelic_pruned_fixed_pruned.prune.in --make-bed --out merged_hg38_biallelic_pruned_fixed_pruned2

# filter individuals 
plink --bfile merged_hg38_biallelic_pruned_fixed_pruned2 --keep keep_plink.txt --make-bed --out merged_hg38_biallelic_pruned_fixed_pruned2_filtered

# you can do this with plink filter. again, I have only done this with VCFs before
# I want to keep all individuals in this list: 
gtex_admixed=/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt
# as well as certain 1000 Genomes populations. I/we need to do multiple iterations with different reference panels, but include the following to start:
# Asian (N = 207)
#   CHB (Han Chinese in Beijing, China)
#   JPT (Japanese in Tokyo, Japan)
# African (N = 405)
#   YRI (Yoruba in Ibadan, Nigeria)
#   GWD (Gambian in Western Divisions in the Gambia)
#   MSL (Mende in Sierra Leone)
#   ESN (Esan in Nigeria)
# European (N = 99)
#   CEU
# so use the subpopulations specified above, but their labels in the .pop file (see below) should be the corresponding superpop (i.e. EUR, AFR, ASN)
# use plink filter in combination with a file that includes subject names of all aforementioned individuals 

# make .pop file to be able to run admixture supervised 
# rows need to be in the same order as the .fam file 
# make one file with SUBJID, SUPERPOP, SUBPOP columns. use '-' as the SUBPOP and SUPERPOP values for GTEx individuals 
# once you make sure SUBJID is in the same order as subjects in the .fam file, extract the SUPERPOP column from that file and give it the same prefix as the bed/bim/fam files with suffix ".pop"

# add software to path 
export PATH=$PATH:/mnt/lab_data/montgomery/nicolerg/TOOLS/admixture_linux-1.3.0
# run admixture 
K=3 # number of reference populations included in the .pop file 
admixture -j24 ${merged_pruned_filtered}.bed ${K} --supervised

paste *pops.txt *.3.Q | sed "s/ /	/g" > ~/gtex-admix/admixture/merged_hg38_biallelic_pruned_fixed_pruned2_filtered.3.Q.subj
# sed -i "s/ /	/g" merged-hg38-thinned-filtered-rfmix-pruned.3.Q.subj



### OLD 

# # LD-prune each merged file 

# gtex_vcf=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz

# # # thin GTEx
# # vcftools --gzvcf ${gtex_vcf} --thin 5000 --recode --out ${outdir}/gtex-hg38-maf01-phased-thinned &

# # # thin 1000 Genomes
# # for chr in {1..22}; do 
# # 	vcftools --gzvcf /mnt/lab_data/montgomery/shared/1KG/hg38/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz --thin 5000 --recode --out ${outdir}/tg-hg38-phased-thinned-chr${chr}
# # done

# # concatenate 1000 Genomes
# vcf-concat ${outdir}/tg-hg38-phased-thinned-chr{1..22}.recode.vcf | bgzip -c > ${outdir}/tg-hg38-phased-thinned.vcf.gz

# wait
# bgzip ${outdir}/gtex-hg38-maf01-phased-thinned.recode.vcf

# # merge GTEx and 1000 Genomes 
# bcftools merge -O z -o ${outdir}/merged-hg38-thinned.vcf.gz --threads 8 ${outdir}/gtex-hg38-maf01-phased-thinned.recode.vcf.gz ${outdir}/tg-hg38-phased-thinned.recode.vcf.gz 

# # filter subjects (make .pop file manually)
# cut -f1 ~/gtex-admix/metadata/admixture.pop.txt > ~/gtex-admix/metadata/keep_admixture.txt # only 100 individuals per superpop (equal random); note some of these are not actually in the VCF
# vcftools --gzvcf ${outdir}/tg-hg38-phased-thinned.recode.vcf.gz  --recode --keep ${outdir}/keep_subj_rfmix.txt --out ${outdir}/tg-hg38-phased-thinned-filtered-rfmix
# #vcftools --gzvcf ${outdir}/tg-hg38-phased-thinned.vcf.gz --recode --keep ~/gtex-admix/metadata/keep_admixture.txt --out ${outdir}/tg-hg38-phased-thinned-filtered

# # tabix -p vcf ${outdir}/gtex-hg38-maf01-phased-thinned-filtered.recode.vcf.gz
# # tabix -p vcf ${outdir}/tg-hg38-phased-thinned-filtered.recode.vcf.gz

# bgzip ${outdir}/tg-hg38-phased-thinned-filtered-rfmix.recode.vcf
# tabix -p vcf ${outdir}/tg-hg38-phased-thinned-filtered-rfmix.recode.vcf.gz

# # now this should be ready to convert to .bed .bim. fam 
# plink --vcf ${outdir}/merged-hg38-thinned-filtered.vcf.gz --make-bed --out ${outdir}/merged-hg38-thinned-filtered

# # order pop file 
# #zgrep -E -m1 "^#CHROM" merged-hg38-thinned-filtered.vcf.gz | sed "s/	/\n/g" > ~/gtex-admix/metadata/vcf_head
# Rscript ~/gtex-admix/scripts/order_pop_file.R 
# #cut -f2 ~/gtex-admix/metadata/merged-hg38-thinned-filtered-pruned.pop > ${outdir}/merged-hg38-thinned-filtered-pruned.pop
# cut -f3 merged-hg38-thinned-filtered-rfmix-pruned-pops.txt > merged-hg38-thinned-filtered-rfmix-pruned.pop

# # run admixture 
# admixture -j12 ${outdir}/merged-hg38-thinned-filtered-pruned.bed 5 --supervised

# paste ../../thinned/merged-hg38-thinned-filtered-rfmix-pruned-pops.txt merged-hg38-thinned-filtered-rfmix-pruned.3.Q > merged-hg38-thinned-filtered-rfmix-pruned.3.Q.subj
# sed -i "s/ /	/g" merged-hg38-thinned-filtered-rfmix-pruned.3.Q.subj



