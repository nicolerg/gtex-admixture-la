# LAVA analysis pipeline


# DEPRECATED - see eqtl-analysis-pipeline 

chr=1
tissue="lung"
indir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/${tissue}/chr${chr}"
scripts="/users/nicolerg/gtex-admix/scripts/local_eqtl"
outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/${tissue}"
picdir="/users/nicolerg/gtex-admix/${tissue}"
allpairs="/mnt/lab_data/montgomery/nicolerg/local-eqtl/${tissue}/Lung.allpairs.chr${chr}.txt.gz"

# collapse chunked files

head -1 ${indir}/${tissue}-local-unkcnt-chr${chr}-00.tsv > ${indir}/${tissue}-local-unkcnt-chr${chr}.tsv
for file in `ls ${indir} | grep "${tissue}-local-unkcnt-chr${chr}-"`
do 
	sed '1d' ${indir}/${file} >> ${indir}/${tissue}-local-unkcnt-chr${chr}.tsv
done

#rm ${indir}/${tissue}-local-unkcnt-chr${chr}-*
gzip "${indir}/${tissue}-local-unkcnt-chr${chr}.tsv"

head -1 ${indir}/${tissue}-local-allpairs-chr${chr}-00.tsv > ${indir}/${tissue}-local-allpairs-chr${chr}.tsv
for file in `ls ${indir} | grep "${tissue}-local-allpairs-chr${chr}-"`
do 
	sed '1d' ${indir}/${file} >> ${indir}/${tissue}-local-allpairs-chr${chr}.tsv
done

#rm ${indir}/${tissue}-local-allpairs-chr${chr}-*
gzip "${indir}/${tissue}-local-allpairs-chr${chr}.tsv"

# merge p-values from LAVA and GTEx

mypairs="${outdir}/lung-local-allpairs-chr${chr}.tsv.gz"
global_local="${outdir}/${tissue}-chr${chr}-gtex-local.tsv"

python ${scripts}/compare-local-results.py ${mypairs} ${allpairs} ${global_local}

gzip ${global_local}
global_local="${global_local}.gz"

# make scatterplot and select outliers

scatter="${picdir}/local-v-gtex-${tissue}-chr${chr}.jpg"
sigpairs="${picdir}/${tissue}-chr${chr}-sigdiff-pairs.tsv"

Rscript ${scripts}/plot-local-results.R ${global_local} ${scatter} ${sigpairs}

# parse sigpairs file

outlocal="${outdir}/${tissue}-chr${chr}-local-sigdifflabels.tsv"
outgtex="${outdir}/${tissue}-chr${chr}-gtex-sigdifflabels.tsv"

# manual curation of mypairs file to remove incomplete lines
mypairs="${outdir}/${tissue}-local-allpairs-chr${chr}-filt.tsv.gz"

python ${scripts}/parse-sigdiff.py ${local_file} ${allpairs} ${sigpairs} ${outlocal} ${outgtex}

# make manhattan plots

Rscript ${scripts}/manhattan.R ${picdir} ${mypairs} ${outlocal} ${outgtex}
