# split expression files first by chromosome and then into chunks of 20 genes

# tissue="wholeblood"
# expm="/mnt/lab_data/montgomery/nicolerg/local-eqtl/${tissue}/Whole_Blood_expression_chr21.bed"
# outdir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/${tissue}"

tissue='lung'
expm='/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/anova/filtered/chr1-expr-filt-byafrsig.tsv.gz'
outdir='/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/chr1/anova/filtered'

#for chr in {1..21}
for chr in 1
do
	zgrep -E "^chr${chr}	" ${expm} > "${outdir}/${tissue}.chr${chr}.expression.bed"
done

#for chr in {1..21}
for chr in 1
do 
	exp="${outdir}/${tissue}.chr${chr}.expression.bed"
	split -d -l 10 ${exp} "${outdir}/chunks/${tissue}.chr${chr}.expression."
done

#zcat ${expm} | head -1 > ${outdir}/header.txt
header='/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/header.txt'
cp ${header} ${outdir}/header.txt

for file in `ls ${outdir}/chunks | grep "expression"`
do
	cat ${header} > tmp.txt
	cat ${outdir}/chunks/${file} >> tmp.txt
	rm ${outdir}/chunks/${file}
	mv tmp.txt ${outdir}/chunks/${file}
done


