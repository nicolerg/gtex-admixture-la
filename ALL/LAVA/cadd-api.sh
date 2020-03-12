#!/bin/bash

tissue=$1
indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}
infile=${indir}/${tissue}.samegene.diffSNP.master.tsv

cut -f3 ${infile} | sed -e '1d' > ${indir}/tmp.txt

while read line; do 
	input=`echo $line | sed "s/_b38//" | sed "s/chr//" | sed "s/_/:/"`
	score=`curl https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.4/${input} | grep "RawScore" | sed "s/\"RawScore\": \"//" | sed "s/\",//"`
	wc=`echo $score | wc -l`
	
	if [ "$wc" -gt 1 ]; then
		continue
	fi
	if [ "$score" == "" ]; then 
		continue
	fi

	echo "${line}	${score}" >> ${indir}/${tissue}.cadd.map.txt
done < ${indir}/tmp.txt

rm ${indir}/tmp.txt

