#!/bin/bash
# run LAVA eQTL calling (LocalAA and GlobalAA)

set -e

module load r/3.5

base=/path/to/analysis/root
srcdir=/path/to/gtex-admixture-la/eqtl
exprdir=${base}/expression_dir # GTEx v8 expression files
localcov=${base}/local_ancestry # local ancestry covariates (e.g. chr10.hg19.maf0.localcov.tsv.gz)
admix_ids=${localcov}/gtex-admixed0.9.txt # see root of repository 
geno=${base}/manuscript_version/admixed_vcf # path to filtered VCFs for eQTL calling 

mkdir -p ${base}/slurm # output folder for log files

for tissue in Whole_Blood; do # list as many tissues as you want; must correspond with GTEx file name prefixes 

    indir=${base}/${tissue}
    edir=${indir}/chunks

    mkdir -p $edir

    exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz
    globalcov=${indir}/${tissue}.all_covariates.txt # path to all covaraites file (GTEx covariates and new PEERs)

    for chr in {1..22}; do
     
        outdir=${indir}/chr${chr}
        mkdir -p ${outdir}
       
        geno=${geno}/gtex.admixed.MAC10.phased.chr${chr}.recode.vcf.gz

        if [[ -f "${indir}/${tissue}.chr${chr}.expression.bed" ]]; then
            echo 'expression already processed'
        else

            # split expression into chunks
            exp=${indir}/${tissue}.chr${chr}.expression.bed
            zgrep -E "^chr${chr}	" ${exprfile} > ${exp}
            split -d -l 50 ${exp} "${edir}/${tissue}.chr${chr}.expression."

            header=${edir}/header.txt
            zcat ${exprfile} | head -1 > ${header}
            for file in `ls ${edir} | grep "chr${chr}.expression"`
            do
                cat ${header} > ${edir}/tmp.chr${chr}.txt
                cat ${edir}/${file} >> ${edir}/tmp.chr${chr}.txt
                rm ${edir}/${file}
                mv ${edir}/tmp.chr${chr}.txt ${edir}/${file}
            done
            gzip ${edir}/${tissue}.chr${chr}.expression*
        fi

        # now run eQTL calling; submit via sbatch

        for chunk in `ls ${edir} | grep "chr${chr}.expression"`; do

            gene=${edir}/${chunk}
            suf=`echo ${gene} | sed 's/.*expression\.//' | sed 's/\..*//'`

            echo "Starting chr $chr chunk $suf"
            out=${outdir}/${tissue}-local-global-allpairs-chr${chr}-${suf}.tsv

            sbatch --export=tissue=${tissue},exprfile=${gene},globalcov=${globalcov},geno=${geno},localcov=${localcov}/chr${chr}.hg19.maf0.localcov.tsv.gz,out=${out},chr=${chr} \
		eqtl_localaa_globalaa_sbatch.sh
        done

    done

done
