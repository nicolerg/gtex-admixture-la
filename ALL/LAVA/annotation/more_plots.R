library(data.table)
library(ggplot2)

finemap = '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/finemap-coloc/filtered_finemap_results.txt'
finemap <- fread(finemap, sep='\t', header=TRUE)

get_method <- function(x){
	if (grepl('_global_',x)){
		return('global')
	} 
	if (grepl('_local_',x)){
		return('local')
	}
}

head(finemap)
finemap[, method := sapply(eqtl_file, get_method)]
dt <- finemap[,.(method, clpp)]

head(dt)

t.test(dt[method=='global',clpp], dt[method=='local',clpp])
ks.test(dt[method=='global',clpp], dt[method=='local',clpp])

jpeg('~/gtex-admix/plots/final_figures/finemap-coloc-clpp.jpg', width=1200, height=600, res=150)
ggplot(dt, aes(fill=factor(method))) +
	geom_histogram(aes(x=as.numeric(clpp)), alpha=0.5, bins=500) +
	scale_fill_manual(values=c(global='blue',local='red'),name='Method') +
	labs(x='CLPP score',y='Frequency') +
	coord_cartesian(xlim=c(0,0.3), ylim=c(0,175))
dev.off()


