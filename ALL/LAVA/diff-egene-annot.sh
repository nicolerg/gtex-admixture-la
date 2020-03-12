#!/bin/bash

# Now I want to look at genes where an eQTL was called in ONLY one method. set significance cutoff

set -e

srcdir=/users/nicolerg/gtex-admix/scripts/LAVA
base=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed

# take output from extract-egene-lava.py

for tissue in Lung Adipose_Subcutaneous Nerve_Tibial Artery_Tibial Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal; do

	indir=${base}/${tissue}
	# compare distributions of N top SNPs (box plot and density plot)
	# output files of SNPs in overlapping eGenes
	global=${indir}/gtex.admix.global.egenes.tied.txt.gz
	lava=${indir}/gtex.admix.lava.egenes.tied.txt.gz
	Rscript ${srcdir}/merge-egenes-2.R ${global} ${lava} ${indir} ${tissue} &

done

# annotate 

for tissue in Lung Adipose_Subcutaneous Nerve_Tibial Artery_Tibial Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal; do

	credir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation"
	cre=`ls ${credir} | grep "$tissue"`

	if [ "$cre" != "" ]; then 
		
		# annotate SNPs

		indir=${base}/${tissue}
		infile=${indir}/${tissue}-uniq-eqtls-ties.txt
		python ${srcdir}/encode-reg-enrich-diffegene.py ${infile} ${credir}/${cre} ${indir} ${tissue} &

	fi
done

# plot annotations 

for tissue in Adipose_Subcutaneous Lung Artery_Tibial Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal; do

	infile=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}/${tissue}.diff.egene.annotated.tsv
	allp=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}/${tissue}-LAVA-global-merged-allpairs-filt.tsv.gz
	Rscript ${srcdir}/plot-annot-pval-2.R ${tissue} ${infile} ${allp} &

done 
