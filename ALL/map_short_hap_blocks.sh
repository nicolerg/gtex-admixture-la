# map short haplotype blocks back to the genome - are a lot of short blocks coming from one place?

cd /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/hap_blocks

zcat /users/nicolerg/gtex-admix/master.ancestry.bed.gz | cut -f 1-3 | sed -e '1d' > tmp 
sort -k1,1 -k2,2n tmp > sorted.tmp.bed

# use R to filter down to intervals < 1MB

bedtools merge -i hapblocks_filt_1MB.bed -c 1 -o count > hapblocks_filt_1MB.bed

# plot result

Rscript ~/gtex-admix/scripts/plot_short_haplo_cov.R
