# regulatory enrichment in genes that are an eQTL with both methods but have different lead SNPs

# for tissue in Lung Muscle_Skeletal; do 

# 	indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}
# 	infile=${indir}/merged.eQTL.diff.leadSNP.tsv
# 	outdir=${indir}

# 	if [ "$tissue" == "Muscle_Skeletal" ]; then
# 		reg=/mnt/lab_data/montgomery/nicolerg/local-eqtl/reg-enrich/homo_sapiens.GRCh38.Psoas_Muscle.Regulatory_Build.regulatory_activity.20161111.gff.gz
# 	elif [ "$tissue" == "Lung" ]; then
# 		reg=/mnt/lab_data/montgomery/nicolerg/local-eqtl/reg-enrich/homo_sapiens.GRCh38.Lung.Regulatory_Build.regulatory_activity.20161111.gff.gz
# 	fi

# 	python ~/gtex-admix/scripts/LAVA/parse-reg-enrich.py ${infile} ${reg} ${outdir} &

# 	wait

# 	global=${indir}/global.merged.eQTL.diff.leadSNP.tsv.gz
# 	lava=${indir}/lava.merged.eQTL.diff.leadSNP.tsv.gz
# 	Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${global} ${lava} ${tissue}

# done

# get p-values from both global and lava for eGenes where lead SNP is different 

for tissue in Lung Muscle_Skeletal; do 

	indir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}
	infile=${indir}/merged.eQTL.diff.leadSNP.tsv
	allpairs=${indir}/${tissue}-LAVA-global-merged-allpairs.tsv.gz
	outfile=${indir}/merged.eQTL.diff.leadSNP.pvalues.tsv

	python ~/gtex-admix/scripts/LAVA/pvalue-pairs.py ${infile} ${allpairs} ${outfile} &

done

