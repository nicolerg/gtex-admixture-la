# plot distribution of short haplotype blocks 

library(data.table)
library(ggplot2)

hap = fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/hap_blocks/hapblocks_filt_1MB_merged_count.txt',sep='\t',header=F)
chr = fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/hap_blocks/chrom_sizes.txt',sep='\t',header=F)
hap = data.table(rbind(hap,chr))

head(hap)
tail(hap)

pdf('~/gtex-admix/plots/august/short-hap-across-genome.pdf',width=20,height=20)
ggplot(hap) + 
	geom_rect(aes(ymin=0,ymax=V4,xmin=as.numeric(V2)/1e6,xmax=as.numeric(V3)/1e6)) +
	facet_wrap(~V1,ncol=4,scales="free_x") + 
	theme_bw() + 
	labs(x='Position (MB)',y='Count of haplotype blocks < 1M bp') 
dev.off()
