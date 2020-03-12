library(ggplot2)

dict <- read.table('/users/nicolerg/gtex-admix/gtexV8-admixture-supp-df.tsv', sep='\t', header=TRUE)
master <- read.table('/users/nicolerg/gtex-admix/gtexV8-rfmix-admixture-df.tsv',sep='\t',header=TRUE)

pdf("/users/nicolerg/gtex-admix/plots/block-all-haplo.pdf", width=8, height=6)
ggplot(master) + geom_segment(aes(x=Start, xend=Stop, y=Y, yend=Y, colour=Pop), size=0.8, lineend="butt") +
	theme_bw() +
	theme(panel.border = element_rect(size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_text(size = 14)) +
	scale_colour_manual(values = c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="gray"), name="Local Ancestry", 
		breaks=c("ASN","EUR","AFR","UNK"),labels=c("Asian","European","African","Unknown")) +
	geom_vline(xintercept=dict$Stop[1:21], size=0.5) + 
	scale_x_discrete(limits=(dict$Start+dict$Stop)/2,
        labels=c(1:22),name="Chromosome") +
  	scale_y_continuous(expand = c(0, 0)) +
  	labs(title="Admixed Individuals (N=114)",y=NULL)
dev.off()

# additional for paper figure

df <- data.frame(c(rep("White",714),rep("African_American",103),rep("Asian",11)))
colnames(df) <- "Race"
df$y <- c(1:nrow(df))
df$xstart <- 0
df$xend <- 1

jpeg("/users/nicolerg/gtex-admix/plots/self-reported-block.jpg", width=200, height=700, res=200)
ggplot(df) + geom_segment(aes(x=xstart, xend=xend, y=y, yend=y, colour=Race), size=0.8, lineend="butt") +
	theme_bw() + 
	theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.position = 'none') +
	scale_colour_manual(values = c(Asian="#ADFF33",White="#0000FF",African_American="#FF9900")) + 
	labs(x = NULL, y = "828 Individuals")
dev.off()
