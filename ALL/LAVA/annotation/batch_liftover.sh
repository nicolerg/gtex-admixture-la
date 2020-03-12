#!/bin/bash 

# liftover hg19 to hg38

chainfile=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode/hg19ToHg38.over.chain.gz

mkdir -p  /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode/old_ref 

while read line; do
	
	hg19=`find "/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode" -name "$line.bed.gz"`
	echo $hg19
	outdir=`dirname $hg19`
	liftOver -bedPlus=4 -tab $hg19 $chainfile ${outdir}/${line}.hg38.bed ${outdir}/${line}.unmapped.bed
	
	mv $hg19 /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode/old_ref/${line}.hg19.bed.gz
	gzip ${outdir}/${line}.hg38.bed ${outdir}/${line}.unmapped.bed
	mv ${outdir}/${line}.unmapped.bed.gz /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode/old_ref/${line}.unmapped.bed.gz

done < /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode/old_ref.txt



# now sort and merge 

base=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode

for tissue in gastrocnemius_medialis subcutaneous_adipose_tissue suprapubic_skin tibial_artery tibial_nerve lung; do 
	indir=${base}/${tissue}
	zcat ${indir}/*.bed.gz | sort -k1,1 -k2,2n | bedtools merge -i stdin > ${indir}/${tissue}.merged.bed 
	gzip ${indir}/${tissue}.merged.bed 
done 

