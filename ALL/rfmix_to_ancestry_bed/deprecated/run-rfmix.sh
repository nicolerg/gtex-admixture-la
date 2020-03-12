# DEPRECATED 

indir="/mnt/lab_data/montgomery/nicolerg/rfmix-admixed"
outdir="/mnt/lab_data/montgomery/nicolerg/rfmix-admixed/whole-outputs"

for chr in {1..22}; do 
	for ind in {81..114}; do
		start=$((${ind}*2-1))
		end=$((${ind}*2))
		cut -d' ' -f${start}-${end},229-1650 ${indir}/classes > ${indir}/classes-${ind} 
		prefix="/mnt/lab_data/montgomery/nicolerg/rfmix-inputs/whole-output/GTEX-${ind}-chr${chr}"
		echo "Working on chromsome ${chr} individual ${ind}..."

		alleles_in=${indir}/alleles-chr${chr}
		cut -c${start}-${end},229-1650 ${alleles_in} > ${indir}/alleles-${ind}-chr${chr}
		
		classes="${indir}/classes-${ind}"
		snp_locations="${indir}/snp_locations-chr${chr}"
		alleles="${indir}/alleles-${ind}-chr${chr}"

		python /users/nicolerg/RFMix_v1.5.4/RunRFMix.py PopPhased ${alleles} ${classes} ${snp_locations} --forward-backward -o ${prefix} &
	done
done
