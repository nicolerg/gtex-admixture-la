set -e 

base=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed

for tissue in Skin_Not_Sun_Exposed_Suprapubic; do

	indir=${base}/${tissue}
	
	for chr in {1..22}; do
		outdir=${indir}/chr${chr}
		f=`ls ${outdir} | grep "${tissue}-LAVA-global-allpairs" | head -1`
		echo "$f"
		head -1 ${outdir}/${f} > ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv
		for chunk in `ls ${outdir} | grep "${tissue}-LAVA-global-allpairs-chr${chr}"`; do
			file=${outdir}/${chunk}
			cat ${file} | sed -e '1d' >> ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv
			rm ${file}
		done
		gzip ${outdir}/${tissue}-LAVA-global-merged-allpairs-chr${chr}.tsv
	done

	outdir=/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/${tissue}

	zcat ${outdir}/chr1/${tissue}-LAVA-global-merged-allpairs-chr1.tsv.gz | head -1 > ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv
	for dir in chr{1..22}; do
		file=${outdir}/${dir}/${tissue}-LAVA-global-merged-allpairs-${dir}.tsv.gz
		zcat $file | sed -e '1d' >> ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv
	done
	gzip ${outdir}/${tissue}-LAVA-global-merged-allpairs.tsv

done
