# Author: Nicole Gay
# Winter 2018
# DEPRECATED
# Note that this script only generates a random (unordered) subset of each chromosome for GTEx and 1000 Genomes separately
# This does not generate subsets of the intersection between GTEx and 1000 Genomes
# Use of this script is NOT recommended to generate PCA input

indir="/mnt/lab_data/montgomery/nicolerg/rfmix-inputs"
outdir="/mnt/lab_data/montgomery/nicolerg/pca-subset"

# first do it for 1000 Genomes

for c in {4..22}
do

	infile="${tgdir}/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
	outfile="${outdir}/1000genomes-chr${c}-subset.vcf"
	zcat ${infile} | perl -ne 'print if (rand() < .001)' > ${outfile}

	sed -i '/^#/d' ${outfile}
	gzip ${outfile}

done

for i in {1..22}; do zcat ${outdir}/1000genomes-chr${i}-subset.vcf.gz >> ${outdir}/1000genomes-rand-all.vcf; done
# prepend "chr"
sed -i "s/^/chr"
gzip ${outdir}/1000genomes-rand-all.vcf

# whoops, need to prepend header
zcat ${tgdir}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -m1 "^#CHROM" > ${outdir}/1000genomes.header.txt
zcat ${outdir}/1000genomes-rand-all.vcf.gz >> ${outdir}/1000genomes.header.txt
mv ${outdir}/1000genomes.header.txt ${outdir}/1000genomes.rand.all.vcf
gzip ${outdir}/1000genomes.rand.all.vcf
rm ${outdir}/1000genomes-rand-all.vcf.gz


# now for GTEx

infile="/mnt/lab_data/montgomery/nicolerg/gtex-hg19/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_hg19.vcf.gz"
outfile="/mnt/lab_data/montgomery/nicolerg/pca-subset/gtex-rand-all.vcf"
zcat ${infile} | perl -ne 'print if (rand() < .002)' > ${outfile}
sed -i '/^#/d' ${outfile}

# need to prepend header
zcat ${infile} | grep -m1 "^#CHROM" > ${outdir}/gtex.header.txt
cat ${outfile} >> ${outdir}/gtex.header.txt
mv ${outdir}/gtex.header.txt ${outdir}/gtex.rand.all.vcf
gzip ${outdir}/gtex.rand.all.vcf
rm ${outfile}