#!/bin/bash

set -e

srcdir=~/gtex-admix/scripts/LAVA

for tissue in Lung Muscle_Skeletal; do

	indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}
	allpairs=${indir}/${tissue}-LAVA-global-merged-allpairs.tsv.gz

	# ======================================================================================

	# # summarize eGenes; define eQTL genes

	# python ${srcdir}/extract-egene-lava.py ${allpairs} ${tissue} ${indir}

	# # # best SNP per gene:
	# # gtex.admix.global.egenes.txt.gz
	# # gtex.admix.lava.egenes.txt.gz

	# ======================================================================================

	# # merge; filter down to genes where lead SNP is different between global and lava

	# outdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}
	# global=${outdir}/gtex.admix.global.egenes.txt.gz
	# lava=${outdir}/gtex.admix.lava.egenes.txt.gz
	# Rscript ${srcdir}/merge-egenes.R ${global} ${lava} ${outdir}

	# ======================================================================================
	 
	# # parse regulatory enrichment 

	# regdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/reg-enrich/

	# if [ "$tissue" == "Muscle_Skeletal" ]; then
	# 	reg=${regdir}/homo_sapiens.GRCh38.Psoas_Muscle.Regulatory_Build.regulatory_activity.20161111.gff.gz
	# elif [ "$tissue" == "Lung" ]; then
	# 	reg=${regdir}/homo_sapiens.GRCh38.Lung.Regulatory_Build.regulatory_activity.20161111.gff.gz
	# fi

	# python ${srcdir}/parse-reg-enrich.py ${indir}/merged-lead-snps.txt ${reg} ${indir} &

done

wait

for tissue in Lung Muscle_Skeletal; do 

	# plot regulatory enrichment 

	indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}

	global=${indir}/global.eqtl.regenrich.tsv.gz
	lava=${indir}/lava.eqtl.regenrich.tsv.gz
	Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${global} ${lava} ${tissue}-eqtl 

	global=${indir}/global.leadsnp.regenrich.tsv.gz
	lava=${indir}/lava.leadsnp.regenrich.tsv.gz
	Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${global} ${lava} ${tissue}-leadsnp

done
