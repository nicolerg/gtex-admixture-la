# Show effect of iteratively adding PCs and PEER factors on eQTL discovery (measured by number of pairs with pval_nominal < 0.05)
# Performed using the nominal pass of FastQTL

base="/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8"
covariates="${base}/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
expression="${base}/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz"
vcf="${base}/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz"
outdir="/mnt/lab_data/montgomery/nicolerg/gtex-eqtl/fastqtl-peer-pc/wholeblood"

###############################################

# source activate py35

# # see how 0 genotype PCs influences discovery

# sed -n -e 1p -e 7,69p ${covariates} > ${outdir}/base.covariates.pc.txt # everything but PCs
# #sed -n -e 1p -e 67,69p ${covariates} > ${outdir}/base.covariates.pc.txt # no PEER or PCs

# taskset -c 0-25 python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py \
# ${vcf} \
# ${expression} \
# "nominal-fastQTL-wholeblood-0pcs" \
# --covariates ${outdir}/base.covariates.pc.txt \
# --window 1e6 \
# --ma_sample_threshold 10 \
# --maf_threshold 0.01 \
# --chunks 100 \
# --threads 36 \
# -o ${outdir} 

# # see how 0 PEER factors influences discovery

# sed -n -e 1,6p -e 67,69p ${covariates} > ${outdir}/base.covariates.peer.txt # everything but PEERs
# #sed -n -e 1p -e 67,69p ${covariates} > ${outdir}/base.covariates.peer.txt # no PEER or PCs

# taskset -c 0-25 python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py \
# ${vcf} \
# ${expression} \
# "nominal-fastQTL-wholeblood-0peer" \
# --covariates ${outdir}/base.covariates.peer.txt \
# --window 1e6 \
# --ma_sample_threshold 10 \
# --maf_threshold 0.01 \
# --chunks 100 \
# --threads 36 \
# -o ${outdir} 

# # no PEER factors or PCs

# outdir="/mnt/lab_data/montgomery/nicolerg/gtex-eqtl/fastqtl-peer-pc/v2"
# sed -n -e 1p -e 67,69p ${covariates} > ${outdir}/base.covariates.txt # no PEER or PCs

# taskset -c 0-25 python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py \
# ${vcf} \
# ${expression} \
# "nominal-fastQTL-wholeblood-0pc0peer" \
# --covariates ${outdir}/base.covariates.txt \
# --window 1e6 \
# --ma_sample_threshold 10 \
# --maf_threshold 0.01 \
# --chunks 100 \
# --threads 36 \
# -o ${outdir} 

# source deactivate

# outdir="/mnt/lab_data/montgomery/nicolerg/gtex-eqtl/fastqtl-peer-pc"
# file1="${outdir}/nominal-fastQTL-wholeblood-0pcs.allpairs.txt.gz"
# file2="${outdir}/nominal-fastQTL-wholeblood-0peer.allpairs.txt.gz"

# zcat ${file1} | awk '{if($4<=0.00001){print $1}}' | sort | uniq | wc -l > ${outdir}/pc-sigcounts-bygene-new.txt
# cat ${outdir}/pc-sigcounts-bygene.txt >> ${outdir}/pc-sigcounts-bygene-new.txt

# zcat ${file2} | awk '{if($4<=0.00001){print $1}}' | sort | uniq | wc -l > ${outdir}/peer-sigcounts-bygene-new.txt
# cat ${outdir}/peer-sigcounts-bygene.txt >> ${outdir}/peer-sigcounts-bygene-new.txt

# outdir="/mnt/lab_data/montgomery/nicolerg/gtex-eqtl/fastqtl-peer-pc/v2"
# file3="${outdir}/nominal-fastQTL-wholeblood-0pc0peer.allpairs.txt.gz"

# zcat ${file3} | awk '{if($4<=0.00001){print $1}}' | sort | uniq | wc -l 

###############################################

source activate py35

# # PCs

# # sequentially add in genotype PCs and see how it influences discovery 
# sed -n -e 1p -e 7,69p ${covariates} > ${outdir}/base.covariates.pc.txt # everything but PCs
# #sed -n -e 1p -e 67,69p ${covariates} > ${outdir}/base.covariates.pc.txt # no PEER or PCs
# sed -n 2,6p ${covariates} > ${outdir}/pc.covs.txt

# # first with no PCs
# taskset -c 0-25 python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py \
# ${vcf} \
# ${expression} \
# "nominal-fastQTL-adipose-0pcs" \
# --covariates ${outdir}/base.covariates.pc.txt \
# --window 1e6 \
# --ma_sample_threshold 10 \
# --maf_threshold 0.01 \
# --chunks 100 \
# --threads 24 \
# -o ${outdir} 

# # now sequentially add a PC
# for i in {1..5}
# do
# 	sed -n ${i}p ${outdir}/pc.covs.txt >> ${outdir}/base.covariates.pc.txt # add a PC
# 	cp ${outdir}/base.covariates.pc.txt ${outdir}/base.covariates.pc.${i}.txt
# 	taskset -c 0-25 python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py \
# 	${vcf} \
# 	${expression} \
# 	"nominal-fastQTL-adipose-${i}pcs" \
# 	--covariates ${outdir}/base.covariates.pc.${i}.txt \
# 	--window 1e6 \
# 	--ma_sample_threshold 10 \
# 	--maf_threshold 0.01 \
# 	--chunks 100 \
# 	--threads 24 \
# 	-o ${outdir} 
# done 

# PEER factors

# sequentially add in PEER factors and see how if affects discovery (# of significant eGenes found)
sed -n -e 1,6p -e 67,69p ${covariates} > ${outdir}/base.covariates.peer.txt # everything but PCs
#sed -n -e 1p -e 67,69p ${covariates} > ${outdir}/base.covariates.peer.txt # no PEER or PCs
sed -n 7,66p ${covariates} > ${outdir}/peer.covs.txt

# first with no PEER factors
taskset -c 24-35 python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py \
${vcf} \
${expression} \
"nominal-fastQTL-adipose-0peer" \
--covariates ${outdir}/base.covariates.peer.txt \
--window 1e6 \
--ma_sample_threshold 10 \
--maf_threshold 0.01 \
--chunks 100 \
--threads 12 \
-o ${outdir} 

for i in {1..60}
do
	sed -n ${i}p ${outdir}/peer.covs.txt >> ${outdir}/base.covariates.peer.txt # add a PEER
	cp ${outdir}/base.covariates.peer.txt ${outdir}/base.covariates.peer.${i}.txt
	taskset -c 24-35 python /users/nicolerg/fastqtl/python/run_FastQTL_threaded.py \
	${vcf} \
	${expression} \
	"nominal-fastQTL-adipose-${i}peer" \
	--covariates ${outdir}/base.covariates.peer.${i}.txt \
	--window 1e6 \
	--ma_sample_threshold 10 \
	--maf_threshold 0.01 \
	--chunks 100 \
	--threads 12 \
	-o ${outdir} 
done

source deactivate

# for file in ${outdir}/nominal-fastQTL-adipose-{0..60}peer.allpairs.txt.gz
# do
# 	zcat ${file} | awk '{if($4<=0.00001){print $1}}' | sort | uniq | wc -l >> ${outdir}/peer-sigcounts-bygene.txt
# done

# for file in ${outdir}/nominal-fastQTL-adipose-{0..5}pcs.allpairs.txt.gz
# do
# 	zcat ${file} | awk '{if($4<=0.00001){print $1}}' | sort | uniq | wc -l >> ${outdir}/pc-sigcounts-bygene.txt
# done
