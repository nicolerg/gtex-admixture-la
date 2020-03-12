library(data.table)
library(ggplot2)

# get distribution of haplotype lengths in ancestry bed file 
# current paths point to a version of RFMix where only 10 individuals were processed
# the purpose of this was to see if using fewer SNPs results in significantly shorter haplotype blocks. it does not 

indir='/mnt/lab_data/montgomery/nicolerg/gtex-tg-merged/rfmix_maf10_20190814/rfmix_outputs/bed/old'

all_lengths <- c()
for (chr in 1:22){
	dt <- fread(sprintf('%s/ancestry.chr%s.bed',indir,chr),sep='\t',header=F)
	colnames(dt) <- c('#chrom','start','end','ancestry','subject_haplotype')

	dt <- dt[order(subject_haplotype, start, decreasing=F)]
	dt[,lengths := end - start]
	dt[,x := rleid(ancestry)]
	lengths <- dt[,list(concat_lengths=sum(lengths)), by='x']
	lengths <- lengths[,concat_lengths]
	all_lengths <- c(lengths, all_lengths)
}
dt <- data.table(x=all_lengths)

pdf('~/gtex-admix/plots/karyo/haplo_lengths_old-1.pdf',width=8,height=4)
ggplot(dt, aes(x/1e6))+
	geom_histogram(bins=500)+
	labs(x='Haplotype block length (MB)',y='Count') +
	xlim(c(0,100))
dev.off()
