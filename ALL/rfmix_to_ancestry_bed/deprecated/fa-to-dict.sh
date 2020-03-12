#!/bin/bash
# Make .dict files for hg38 chromosomes
for n in {1..22}
do 
	input="/mnt/lab_data/montgomery/nicolerg/hg38-ref/chr${n}.fa"
	output="/mnt/lab_data/montgomery/nicolerg/hg38-ref/chr${n}.dict"
	java -jar picard.jar CreateSequenceDictionary R=${input} O=${output}
done