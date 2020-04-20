#!/bin/bash
# run LAVA eQTL calling (LocalAA and GlobalAA) with permutations to account for multiple testing
# expression is already processed from eQTL calling 

set -e

module load r/3.5

base=/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE
srcdir=/oak/stanford/groups/smontgom/nicolerg/src/gtex-admixture-la/eqtl
#exprdir=/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices # GTEx v8 expression files
localcov=/oak/stanford/groups/smontgom/nicolerg/LAVA/manuscript_version/admixed_ancestry # local ancestry covariates (e.g. chr10.hg19.maf0.localcov.tsv.gz)
admix_ids=/oak/stanford/groups/smontgom/nicolerg/src/gtex-admixture-la/metadata/gtex-admixed0.9.txt 
geno_dir=/oak/stanford/groups/smontgom/nicolerg/LAVA/manuscript_version/admixed_vcf # path to filtered VCFs for eQTL calling 
egenes_master=/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISED_PERMUTE/egenes_master.RData 

mkdir -p ${base}/slurm # output folder for log files

for tissue in Muscle_Skeletal Adipose_Subcutaneous Lung; do # list as many tissues as you want; must correspond with GTEx file name prefixes 

    # indir=${base}/${tissue}
    # edir=${indir}/chunks
    edir=/oak/stanford/groups/smontgom/nicolerg/LAVA/manuscript_version/eqtl_calling/${tissue}/chunks

    #mkdir -p $edir 

    #exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz
    globalcov=${base}/${tissue}/${tissue}.all_covariates.txt # path to all covariates file (GTEx covariates and new PEERs)

    #for chr in {1..22}; do
    for chr in 22; do 
     
        outdir=${base}/${tissue}/chr${chr}
        # outdir=${indir}/chr${chr}
        mkdir -p ${outdir}
       
        geno=${geno_dir}/gtex.admixed.MAC10.phased.chr${chr}.recode.vcf.gz

        # if [[ -f "${indir}/${tissue}.chr${chr}.expression.bed" ]]; then
        #     echo 'expression already processed'
        # else

        #     # split expression into chunks
        #     exp=${indir}/${tissue}.chr${chr}.expression.bed
        #     zgrep -E "^chr${chr}	" ${exprfile} > ${exp}
        #     split -d -l 50 ${exp} "${edir}/${tissue}.chr${chr}.expression."

        #     header=${edir}/header.txt
        #     zcat ${exprfile} | head -1 > ${header}
        #     for file in `ls ${edir} | grep "chr${chr}.expression"`
        #     do
        #         cat ${header} > ${edir}/tmp.chr${chr}.txt
        #         cat ${edir}/${file} >> ${edir}/tmp.chr${chr}.txt
        #         rm ${edir}/${file}
        #         mv ${edir}/tmp.chr${chr}.txt ${edir}/${file}
        #     done
        #     gzip ${edir}/${tissue}.chr${chr}.expression*
        # fi

        # now run permutation procedure; submit via sbatch

        for chunk in `ls ${edir} | grep "chr${chr}.expression"`; do

            gene=${edir}/${chunk}
            suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`

            echo "Starting chr $chr chunk $suf"
            out=${outdir}/${tissue}-local-global-permuted-chr${chr}-${suf}.tsv

            sbatch --export=tissue=${tissue},exprfile=${gene},globalcov=${globalcov},geno=${geno},localcov=${localcov}/chr${chr}.hg19.maf0.localcov.tsv.gz,out=${out},chr=${chr},egenes_master=${egenes_master},subjid=${admix_ids} /oak/stanford/groups/smontgom/nicolerg/src/gtex-admixture-la/eqtl/slurm/eqtl_localaa_globalaa_permute_sbatch.sh

        done

    done

done
