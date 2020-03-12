library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
indir <- args[1]
tissue_rep <- args[2]
df <- args[3]

indir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/"

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# plot replication of GTEx in GEUVADIS ========================================================================================================================

replication <- read.csv(tissue_rep, header=TRUE, sep='\t', stringsAsFactors=FALSE)

head(replication)

replication$EUR_frac <- as.numeric(replication$EUR_frac)*100
replication$AFR_frac <- as.numeric(replication$AFR_frac)*100
replication <- replication[order(replication$EUR_frac, decreasing=TRUE),]

tissuecols <- c(Adipose_Subcutaneous="#FF6600",
	Adipose_Visceral_Omentum="#FFAA00",
	Adrenal_Gland="#33DD33",
	Artery_Aorta="#FF5555",
	Artery_Coronary="#FFAA99",
	Artery_Tibial="#FF0000",
	Brain_Amygdala="#EEEE00",
	Brain_Anterior_cingulate_cortex_BA24="#EEEE00",
	Brain_Caudate_basal_ganglia="#EEEE00",
	Brain_Cerebellar_Hemisphere="#EEEE00",
	Brain_Cerebellum="#EEEE00",
	Brain_Cortex="#EEEE00",
	Brain_Frontal_Cortex_BA9="#EEEE00",
	Brain_Hippocampus="#EEEE00",
	Brain_Hypothalamus="#EEEE00",
	Brain_Nucleus_accumbens_basal_ganglia="#EEEE00",
	Brain_Putamen_basal_ganglia="#EEEE00",
	"Brain_Spinal_cord_cervical_c-1"="#EEEE00",
	Brain_Substantia_nigra="#EEEE00",
	Breast_Mammary_Tissue="#33CCCC",
	"Cells_EBV-transformed_lymphocytes"="#CC66FF",
	Cells_Cultured_fibroblasts="#AAEEFF",
	Colon_Sigmoid="#EEBB77",
	Colon_Transverse="#CC9955",
	Esophagus_Gastroesophageal_Junction="#8B7355",
	Esophagus_Mucosa="#552200",
	Esophagus_Muscularis="#BB9988",
	Heart_Atrial_Appendage="#9900FF",
	Heart_Left_Ventricle="#660099",
	Liver="#AABB66",
	Lung="#99FF00",
	Kidney_Cortex="yellow",
	Minor_Salivary_Gland="gray",
	Muscle_Skeletal="#AAAAFF",
	Nerve_Tibial="#FFD700",
	Ovary="#FFAAFF",
	Pancreas="#995522",
	Pituitary="#AAFF99",
	Prostate="#DDDDDD",
	Skin_Not_Sun_Exposed_Suprapubic="#0000FF",
	Skin_Sun_Exposed_Lower_leg="#7777FF",
	Small_Intestine_Terminal_Ileum="#555522",
	Spleen="#778855",
	Stomach="#FFDD99",
	Testis="#AAAAAA",
	Thyroid="#006600",
	Uterus="#FF66FF",
	Vagina="#FF5599",
	Whole_Blood="#FF00BB")

replication$label <- ""
replication$label[replication$tissue == "Cells_EBV-transformed_lymphocytes"] <- "Cells_EBV-transformed_lymphocytes"
replication$label[replication$tissue == "Spleen"] <- "Spleen"
replication$label[replication$tissue == "Muscle_Skeletal"] <- "Muscle_Skeletal"
replication$label[replication$tissue == "Small_Intestine_Terminal_Ileum"] <- "Small_Intestine_Terminal_Ileum"

g1 <- ggplot(replication, aes(x=tissue,y=EUR_frac,fill=factor(tissue))) +
	geom_bar(stat="identity") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),axis.title=element_text(size=10), legend.position="none") +
	labs(x="GTEx v8 tissue", y="% of GEUVADIS EUR eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	scale_fill_manual(values=tissuecols) +
	geom_text(aes(label=sample_size),position=position_dodge(width=0.9),vjust=-0.5,size=2) 
g2 <- ggplot(replication, aes(x=sample_size, y=EUR_frac, colour=factor(tissue), label=label)) +
	geom_point() +
	theme_bw() +
	scale_colour_manual(values=tissuecols, name="Tissue") +
	labs(x="Sample size",y="% of GEUVADIS EUR eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	geom_text(aes(label=label),hjust=0, vjust=0) +
	theme(legend.position="none") +
	xlim(c(min(replication$sample_size),800))

pdf("~/gtex-admix/plots/geu-gtex-eur-eqtl-replication-bytissue.pdf", width=8, height=14)
multiplot(g1,g2,cols=1)
dev.off()

replication <- replication[order(replication$AFR_frac, decreasing=TRUE),]

g1 <- ggplot(replication, aes(x=tissue,y=AFR_frac,fill=factor(tissue))) +
	geom_bar(stat="identity") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),axis.title=element_text(size=10), legend.position="none") +
	labs(x="GTEx v8 tissue", y="% of GEUVADIS AFR eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	scale_fill_manual(values=tissuecols) +
	geom_text(aes(label=sample_size),position=position_dodge(width=0.9),vjust=-0.5,size=2) 
g2 <- ggplot(replication, aes(x=sample_size, y=AFR_frac, colour=factor(tissue), label=label)) +
	geom_point() +
	theme_bw() +
	scale_colour_manual(values=tissuecols, name="Tissue") +
	labs(x="Sample size",y="% of GEUVADIS AFR eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	geom_text(aes(label=label),hjust=0, vjust=0) +
	theme(legend.position="none") +
	xlim(c(min(replication$sample_size),800))

pdf("~/gtex-admix/plots/geu-gtex-afr-eqtl-replication-bytissue.pdf", width=8, height=14)
multiplot(g1,g2,cols=1)
dev.off()

# plot replication in LCLs =============================================================

df <- read.delim(df, sep='\t', header=TRUE, stringsAsFactors=FALSE)

df$FRAC <- 100*df$FRAC
df$POP <- paste0(df$POP, " SNPs")

pdf("~/gtex-admix/plots/geu-gtex-lcl-replication.pdf", width=6, height=5)
ggplot(df, aes(x=POP,y=FRAC,group=GROUP,fill=GROUP)) +
	geom_bar(position=position_dodge(),stat="identity") +
	scale_fill_manual(name='Effect Size',values=c("higheur"="blue","highafr"="orange"),labels=c("higheur"="EUR > AFR (N = 32)","highafr"="AFR > EUR (N = 33)")) +
	labs(x="",y="% of GEU eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	theme_bw() +
	theme(axis.text.x = element_text(size=14),axis.title=element_text(size=11)) 
dev.off()

# plot sample size ======================================================================

pheno <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

pheno <- fread(pheno, sep='\t', header=TRUE)
aa <- pheno[RACE==2,SUBJID] # only AA

flist <- dir(path = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices", 
	pattern = "normalized_expression.bed.gz", 
	full.names = TRUE)

t <- c()
n <- c()
for (f in flist){
	if (!grepl("gz.tbi",f,fixed=TRUE)){
		x <- fread(sprintf("zcat %s",f),sep='\t',header=TRUE,nrows=3)
		tissue <- gsub(".*/","",f)
		tissue <- gsub(".v8.normalized_expression.bed.gz","",tissue,fixed=TRUE)
		print(tissue)
		subjs <- colnames(x)
		subjs <- subjs[subjs %in% aa]
		naa <- length(subjs)
		print(naa)
		t <- c(t, tissue)
		n <- c(n, naa)
	}
}

df <- data.frame(tissue=t,N_AA=n)

head(df)

pdf("~/gtex-admix/plots/gtex-AA-pertissue.pdf", width=8, height=6)
ggplot(df,aes(x=tissue,y=N_AA,fill=tissue)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=tissuecols,name="GTEx Tissue") +
	labs(x="GTEx Tissue",y="N African American samples") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),axis.title=element_text(size=10), legend.position="none") +
	geom_hline(yintercept=57,linetype="dashed")
dev.off()

# plot sample size versus replication ======================================================================

replication$label <- ""
replication$label[replication$tissue == "Cells_EBV-transformed_lymphocytes"] <- "Cells_EBV-transformed_lymphocytes"
# replication$label[replication$tissue == "Lung"] <- "Lung"
# replication$label[replication$tissue == "Esophagus_Mucosa"] <- "Esophagus_Mucosa"
# replication$label[replication$tissue == "Nerve_Tibial"] <- "Nerve_Tibial"
# replication$label[replication$tissue == "Skin_Not_Sun_Exposed_Suprapubicn"] <- "Skin_Not_Sun_Exposed_Suprapubic"

scat <- merge(replication, df, by="tissue")

pdf("~/gtex-admix/plots/gtex-AA-by-replication.pdf", width=8, height=6)
ggplot(scat, aes(x=N_AA, y=AFR_frac, colour=factor(tissue), label=label)) +
	geom_point() +
	theme_bw() +
	scale_colour_manual(values=tissuecols) +
	labs(x="No. of AA individuals",y="% of GEUVADIS AFR eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	geom_text(aes(label=label),hjust=0, vjust=0) +
	theme(legend.position="none") +
	geom_vline(xintercept=57,linetype="dashed")
dev.off()


