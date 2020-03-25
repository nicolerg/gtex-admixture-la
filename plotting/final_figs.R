# final figures
library(data.table)
library(ggplot2)
library(gridExtra)
library(reshape)
library(OneR)
library(ggrepel)
library(grid)
library(gridExtra)

source('~/gtex-admix/scripts/LAVA/individual_final_figs.R')

# # Genome Biology figure standards

# Each figure should be closely cropped to minimize the amount of white space surrounding the illustration. 
# Cropping figures improves accuracy when placing the figure in combination with other elements when the accepted manuscript is prepared for publication on our site. 
# width of 600 pixels (standard), 1200 pixels (high resolution).
# width of 85 mm for half page width figure
# width of 170 mm for full page width figure
# maximum height of 225 mm for figure and legend
# image resolution of approximately 300 dpi (dots per inch) at the final size
# Figures should be designed such that all information, including text, is legible at these dimensions. 
# All lines should be wider than 0.25 pt when constrained to standard figure widths. All fonts must be embedded.

# jpeg(width=1200, height=, quality=100, res=300)
# pdf(width=6.69291, height=8.85827)

g_legend <- function(a.gplot){
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}

system("mkdir -p ~/gtex-admix/plots/manuscript")


figure1 <- function(){

	pdf("~/gtex-admix/plots/manuscript/figure1.pdf",width=7.5, height=7.5)
	grid.arrange(genotype_pc(), admix_per_tissue(), rfmix_pc_cor(),
		layout_matrix = rbind(c(1,2),
								c(1,2),
								c(3,2))
	)
	dev.off()
}
figure1()

# add haplo block manually

figure2 <- function(){
	
	#r <- rectGrob(gp=gpar(fill="white",colour='white'))

	# jpeg("~/gtex-admix/plots/manuscript/figure2_1e-05.jpg",width=10,height=6,quality=100,res=300,units='in')
	# grid.arrange(egene_discovery(cutoff=1e-05), qq, venn(cutoff=1e-05), pval_distn(cutoff=1e-05),
	# 	layout_matrix = rbind(c(2,2,2,1,1,1,1),
	# 						c(3,3,3,3,3,4,4))
	# )
	# dev.off()

	pdf("~/gtex-admix/plots/manuscript/figure2_1e-06.pdf",width=7.5,height=6)
	grid.arrange(egene_discovery(cutoff=1e-06), qq(), venn(cutoff=1e-06), pval_distn(cutoff=1e-06),
		layout_matrix = rbind(c(2,2,2,1,1,1,1,1),
							c(2,2,2,1,1,1,1,1),
							c(3,3,3,3,3,4,4,4),
							c(3,3,3,3,3,4,4,4),
							c(3,3,3,3,3,4,4,4))
	)
	dev.off()

}

figure3 <- function(){

	local_var <- local_var_with_fst()
	gtex <- gtex_r2()
	high <- gtex_esnps_r2_high()

	pdf("~/gtex-admix/plots/manuscript/figure3.pdf",width=7.5, height=7)
	grid.arrange(local_var, gtex, high,
		layout_matrix = rbind(c(1,1,1,1,1,1),
							c(1,1,1,1,1,1),
							c(2,2,3,3,3,3),
							c(2,2,3,3,3,3),
							c(2,2,3,3,3,3))
	)
	dev.off()

}

figure4 <- function(cutoff=1e-6){

	coloc <- coloc()
	pop_fst <- fst_distn()
	max_pop <- max_fst_pop()
	r <- rectGrob(gp=gpar(fill="white",col='white'))

	pdf("~/gtex-admix/plots/manuscript/figure4.pdf",width=11, height=12)
	grid.arrange(coloc, r, max_pop, pop_fst, 
		layout_matrix=rbind(c(1,1,1,1,1),
							c(1,1,1,1,1),
							c(2,2,2,3,3),
							c(4,4,4,3,3))
	)
	dev.off()

}


