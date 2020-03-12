srcdir="/users/nicolerg/gtex-admix/scripts/geuvadis"
indir="/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl"
outdir=${indir}

# python ${srcdir}/rsid-hg38-coord.py 

# gtexdir="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations"
# for file in `ls ${gtexdir}`; do 
# 	echo `echo ${file} | sed "s/\..*//"`
# 	python ${srcdir}/geu-gtex-tissue-rep.py ${file} ${outdir} &

# 	running=`ps -ef | grep "nicolerg" | grep "geu-gtex-tissue-rep" | wc -l`
# 	while [ $running -gt 10 ]; do
# 		sleep 180
# 		running=`ps -ef | grep "nicolerg" | grep "geu-gtex-tissue-rep" | wc -l`
# 	done
# done

# wait

# # # merge files
# merged_rep="${outdir}/geu.gtex.rep.merged.txt"
# echo "tissue"	"EUR_frac"	"AFR_frac"	"sample_size" > ${merged_rep}
# for file in `ls ${outdir} | grep "geu.gtex.replication"`; do 
# 	cat $file >> ${outdir}/geu.gtex.rep.merged.txt
# done

# # figure 3 

# Rscript ${srcdir}/geu-effect-sizes.R 
# # GEU-eQTLS-highAFR-effect.tsv
# # GEU-eQTLS-highEUR-effect.tsv

# # figure 4, 1, 2

# python ${srcdir}/geu-eQTL-gtex-lcl-replication.py 

merged_rep="${outdir}/geu.gtex.rep.merged.txt"
df=${indir}/GEU-GTEX-LCL-diff-rvalue-plotdf.tsv
Rscript ${srcdir}/geu-gtex-compare.R ${indir} ${merged_rep} ${df}
# df is from preceding python script
# merged_rep is from geu-gtex-tissue-rep.py

# figure 5

# need to run LAVA and global eQTL calling in a tissue with high eQTL overlap in GEU but also high number of AA individuals




