#!/bin/bash
# run LAVA eQTL calling (LocalAA and GlobalAA)
# for new PEER factors, first run new_peer_admixed.sh to generate covariate files 

set -e

module load r/3.5

base=/labs/smontgom/nicolerg/LAVA
srcdir=/labs/smontgom/nicolerg/src/gtex-admixture-la/eqtl
exprdir=${base}/Whole_Blood
localcov=${base}/manuscript_version/admixed_ancestry
admix_ids=${localcov}/gtex-admixed0.9.txt
geno=${base}/manuscript_version/admixed_vcf

mkdir -p ${base}/slurm

for tissue in Whole_Blood; do 

    indir=${base}/${tissue}
    edir=${indir}/chunks

    mkdir -p $edir

    exprfile=${exprdir}/${tissue}.v8.normalized_expression.bed.gz
    globalcov=${indir}/${tissue}.all_covariates.txt

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
            out=${outdir}/${tissue}-LAVA-global-allpairs-chr${chr}-${suf}.tsv

            sbatch --export=tissue=${tissue},exprfile=${gene},globalcov=${globalcov},geno=${geno},localcov=${localcov}/chr${chr}.hg19.maf0.localcov.tsv.gz,out=${out},chr=${chr} \
		eqtl_localaa_globalaa_sbatch.sh
        done

    done

done
