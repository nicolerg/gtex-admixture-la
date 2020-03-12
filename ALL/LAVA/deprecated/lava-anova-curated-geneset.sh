indir="/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl"
outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/geuvadis-curated"

chr=$1
tissue="Artery_Tibial"

# generate chunks from the curated gene list
geneset="${indir}/eqtl.curated.geneset.chr${chr}.tsv"

split -d -l 5 ${geneset} "${indir}/chunks/eqtl.curated.geneset.chr${chr}."

for file in `ls ${indir}/chunks | grep "eqtl.curated.geneset.chr${chr}"`; do
	head -1 ${geneset} > ${indir}/chunks/tmp.txt
	cat ${indir}/chunks/${file} >> ${indir}/chunks/tmp.txt
	rm ${indir}/chunks/${file}
	mv ${indir}/chunks/tmp.txt ${indir}/chunks/${file}
done

for chunk in `ls $indir/chunks | grep "eqtl.curated.geneset.chr${chr}"`; do
	geneset=${indir}/chunks/${chunk}
	suf=`echo ${geneset} | sed 's/.*\.//' `
	echo "Starting chunk $suf"
	out1="${outdir}/${tissue}-curated-lava-anova-allpairs-chr${chr}-${suf}.tsv"
	taskset -c 30-47 Rscript /users/nicolerg/gtex-admix/scripts/LAVA/lava-anova-curated-geneset-tiss.R ${chr} ${geneset} ${out1} ${tissue} &
	
	running=`ps -ef | grep "nicolerg" | grep "lava-anova-curated-geneset" | wc -l`
	while [ $running -gt 10 ] # run max 15 jobs at a time
	do
		sleep 180
		running=`ps -ef | grep "nicolerg" | grep "lava-anova-curated-geneset" | wc -l`
	done
done

wait

# # merge chunks
# for chunk in `ls $outdir | grep "${tissue}-curated-lava-anova-allpairs-chr${chr}-"`; do
# 	file=${outdir}/${chunk}
# 	cat $file >> ${outdir}/${tissue}-curated-lava-anova-allpairs-chr${chr}.merged.tsv
# done
# rm ${outdir}/${tissue}-curated-lava-anova-allpairs-chr${chr}-*
# gzip ${outdir}/${tissue}-curated-lava-anova-allpairs-chr${chr}.merged.tsv

