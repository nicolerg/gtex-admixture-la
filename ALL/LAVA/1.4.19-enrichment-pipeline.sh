#!/bin/bash

#########################

# ${tissue}-LAVA-global-merged-allpairs.tsv.gz --> extract-egene-lava.py --> gtex.admix.[lava,global].egenes.tied.txt.gz
#	gtex.admix.[lava,global].egenes.tied.txt.gz includes all top SNPs per each TESTED gene (not just eGenes)

# gtex.admix.[lava,global].egenes.tied.txt.gz --> merge-egenes.R --> ${tissue}-[lava,global]-eqtls-ties.txt
#	${tissue}-[lava,global]-eqtls-ties.txt includes ONLY eGenes but ALL SNPs (including same per method)
# gtex.admix.[lava,global].egenes.tied.txt.gz --> encode-reg-enrich.py -->

# gtex.admix.[lava,global].egenes.tied.txt.gz --> encode-reg-enrich.py --> ${tissue}.diff.SNP.annot.tsv, ${tissue}.diff.eSNP.annotated.tsv
# 	${tissue}.diff.SNP.annot.tsv --> gene-level summary of both methods
# 	${tissue}.diff.eSNP.annotated.tsv --> one line per SNP per method per gene

# ignore this one for now: 
# # gtex.admix.[lava,global].egenes.tied.txt.gz --> encode-reg-enrich.py --> ${tissue}.[lava,global].annot.prob.tsv; ${tissue}.[global,lava].eSNP.annotated
# #	${tissue}.[lava,global].annot.prob.tsv = fraction of SNPs in each gene-level top set that have an annotation
# #	${tissue}.[global,lava].eSNP.annotated = annotation of every individual SNP 

# ${tissue}.[lava,global].annot.prob.tsv --> plot-annot-pval.R --> 


set -e

tissue=$1

srcdir=/users/nicolerg/gtex-admix/scripts/LAVA
base=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed
indir=${base}/${tissue}

# =======================

# get top SNPs for each method 
infile=${indir}/${tissue}-LAVA-global-merged-allpairs.tsv.gz
python ${srcdir}/extract-egene-lava.py ${infile} ${tissue} ${indir}

# =======================

for tissue in Lung Adipose_Subcutaneous Nerve_Tibial Artery_Tibial Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal; do

	indir=${base}/${tissue}
	# compare distributions of N top SNPs (box plot and density plot)
	# output files of SNPs in overlapping eGenes
	global=${indir}/gtex.admix.global.egenes.tied.txt.gz
	lava=${indir}/gtex.admix.lava.egenes.tied.txt.gz
	Rscript ${srcdir}/merge-egenes.R ${global} ${lava} ${indir} ${tissue} &

done

# =======================

# annotate SNPs - ENCODE Candidate Regulatory Elements

for tissue in Lung Adipose_Subcutaneous Nerve_Tibial Artery_Tibial Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal; do

	credir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation"
	cre=`ls ${credir} | grep "$tissue"`

	if [ "$cre" != "" ]; then 
		
		# annotate SNPs

		indir=${base}/${tissue}
		infile=${indir}/${tissue}.samegene.diffSNP.master.tsv
		python ${srcdir}/encode-reg-enrich-difflead.py ${infile} ${credir}/${cre} ${indir} ${tissue} &

	fi
done

# filter allpairs files for plotting purposes

for tissue in Lung Adipose_Subcutaneous Nerve_Tibial Artery_Tibial Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal; do

	python ~/gtex-admix/scripts/LAVA/filter_allpairs.py /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}/${tissue}-LAVA-global-merged-allpairs.tsv.gz &

done

# plot annotation 

for tissue in Adipose_Subcutaneous Lung Artery_Tibial Nerve_Tibial Skin_Not_Sun_Exposed_Suprapubic Muscle_Skeletal; do

	infile=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}/${tissue}.diff.SNP.annot.tsv
	full=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}/${tissue}.diff.eSNP.annotated.tsv
	allp=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}/${tissue}-LAVA-global-merged-allpairs-filt.tsv.gz
	caddmap=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}/${tissue}.cadd.map.txt
	Rscript ${srcdir}/plot-annot-pval.R ${infile} ${tissue} ${full} ${allp} ${caddmap} &

done 

