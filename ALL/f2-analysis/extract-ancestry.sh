# DO NOT USE - now prepended to compare-f2-1kg.sh

tmp="/mnt/lab_data/montgomery/nicolerg/tmp"
outdir="/mnt/lab_data/montgomery/nicolerg/rfmix-outputs/f2-variants"

chr="chr22"

rm "${outdir}/bedintersect.${chr}.out"

echo "#CHROM	START	SNP_ID	REF	ALT	SUBJ1	SUBJ2	POP" > "${outdir}/bedintersect.${chr}.out"

counter=1
while read p
do
	echo $p > "${tmp}/anc.${chr}.line.${counter}.bed"
	sed -i 's/ /	/g' "${tmp}/anc.${chr}.line.${counter}.bed"

	ancbed="/mnt/lab_data/montgomery/nicolerg/rfmix-outputs/haplotype-bed/ancestry.${chr}.bed"
	# subjects IDs in fields 7 and 8
	subj1=`cut -f7 "${tmp}/anc.${chr}.line.${counter}.bed"`
	subj2=`cut -f8 "${tmp}/anc.${chr}.line.${counter}.bed"`
	# filter BED file
	grep -E "${subj1}|${subj2}" ${ancbed} > "${tmp}/anc.${chr}.filt.${counter}.bed"

	# run bedtools intersect
	bedtools intersect -a "${tmp}/anc.${chr}.line.${counter}.bed" -b "${tmp}/anc.${chr}.filt.${counter}.bed" -wo > "${tmp}/anc.${chr}.bedintersect"

	# check if line count is at least 2
	if [ `wc -l "${tmp}/anc.${chr}.bedintersect" | sed 's/ .*//'` -gt 1 ]; then
		# check if populations are the same
		pop1=`cut -f12 "${tmp}/anc.${chr}.bedintersect" | head -1`
		pop2=`cut -f12 "${tmp}/anc.${chr}.bedintersect" | tail -1`
		if [ "$pop1" = "$pop2" ]; then
			head -1 "${tmp}/anc.${chr}.bedintersect" | cut -f 1,2,4,5,6,7,8,12 >> "${outdir}/bedintersect.${chr}.out"
		fi
	fi

	rm "${tmp}/anc.${chr}.line.${counter}.bed" "${tmp}/anc.${chr}.filt.${counter}.bed" "${tmp}/anc.${chr}.bedintersect"

done < "/mnt/lab_data/montgomery/nicolerg/rfmix-outputs/f2-variants/sorted.chr22.f2.var.bed"

sed -i "s/^chr//" "${outdir}/bedintersect.${chr}.out"
