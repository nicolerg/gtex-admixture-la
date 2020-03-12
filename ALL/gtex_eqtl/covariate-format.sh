# Author: Nicole Gay
# Winter 2018
# Add sex, platform, and pcr to covariate files
# Reformat covariate files with cov-format.R

outdir="/mnt/lab_data/montgomery/nicolerg/eqtl"
output="/mnt/lab_data/montgomery/nicolerg/eqtl/other_covariates.txt"
temp1=${outdir}/temp1.txt
temp2=${outdir}/temp2.txt
temp3=${outdir}/temp3.txt

cov="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates"

for file in `ls ${cov}`
do 
	if [ ! -e ${temp1} ]
	then
		echo "first sample"
		head -1 ${cov}/${file} > ${temp1}
		tail -3 ${cov}/${file} >> ${temp1}
	else
		echo ${file}
		head -1 ${cov}/${file} > ${temp2}
		tail -3 ${cov}/${file} >> ${temp2}
		paste ${temp1} ${temp2} > ${temp3}
		rm ${temp1}
		mv ${temp3} ${temp1}
	fi
done

mv ${temp1} ${output}
rm ${temp2}

for file in `ls ${outdir} | grep "PEER_covariates"`
do
    prefix=`echo $file | sed 's/\..*//'`
    outfile="${outdir}/${prefix}.all_covariates.txt"
    Rscript concat-covariates.R ${output} ${outdir}/${file} ${outfile} 
done
