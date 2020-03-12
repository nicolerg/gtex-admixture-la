library(data.table)
library(ggplot2)

indir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/"

# plot replication of GTEx in GEUVADIS ========================================================================================================================

replication <- read.csv(paste0(indir,'/mnt/lab_data/montgomery/nicolerg/geuvadis/published-eqtl/'), row.names=1, header=TRUE, sep='\t', stringsAsFactors=FALSE)
replication$tissue <- rownames(replication)
rownames(replication) <- NULL

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
	Breast_Mammary_Tissue="#33CCCC",
	"Cells_EBV-transformed_lymphocytes"="#CC66FF",
	Cells_Transformed_fibroblasts="#AAEEFF",
	Colon_Sigmoid="#EEBB77",
	Colon_Transverse="#CC9955",
	Esophagus_Gastroesophageal_Junction="#8B7355",
	Esophagus_Mucosa="#552200",
	Esophagus_Muscularis="#BB9988",
	Heart_Atrial_Appendage="#9900FF",
	Heart_Left_Ventricle="#660099",
	Liver="#AABB66",
	Lung="#99FF00",
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

pdf("~/gtex-admix/plots/geu-gtex-eur-eqtl-replication-bytissue.pdf", width=12, height=6)
ggplot(replication, aes(x=tissue,y=EUR_frac)) +
	geom_bar(stat="identity",fill=factor(tissue)) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title=element_text(size=12), legend.position="none") +
	theme_bw() +
	labs(x="GTEx v8 tissue", y="% of GEUVADIS EUR eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	scale_fill_manual(values=tissuecols)
ggplot(replication, aes(x=samp_size, y=EUR_frac, colour=factor(tissue))) +
	geom_point() +
	theme_bw() +
	scale_fill_manual(values=tissuecols, name="Tissue") +
	labs(x="Sample size",y="% of GEUVADIS EUR eQTLs with GTEx SNP-gene pairs (p < 0.05)")
dev.off()

replication <- replication[order(replication$EUR_frac, decreasing=TRUE),]

pdf("~/gtex-admix/plots/geu-gtex-afr-eqtl-replication-bytissue.pdf", width=12, height=6)
ggplot(replication, aes(x=tissue,y=AFR_frac)) +
	geom_bar(stat="identity",fill=factor(tissue)) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title=element_text(size=12), legend.position="none") +
	theme_bw() +
	labs(x="GTEx v8 tissue", y="% of GEUVADIS AFR eQTLs with GTEx SNP-gene pairs (p < 0.05)") +
	scale_fill_manual(values=tissuecols)
ggplot(replication, aes(x=samp_size, y=AFR_frac, colour=factor(tissue))) +
	geom_point() +
	theme_bw() +
	scale_fill_manual(values=tissuecols, name="Tissue") +
	labs(x="Sample size",y="% of GEUVADIS AFR eQTLs with GTEx SNP-gene pairs (p < 0.05)")
dev.off()

# plot effect size in GEUVADIS eQTLs ========================================================================================================================

eurgeu <- "EUR373.gene.cis.FDR5.best.rs137.txt.gz"
afrgeu <- "YRI89.gene.cis.FDR5.best.rs137.txt.gz"

eur <- fread(input=sprintf("zcat %s",eurgeu), sep='\t', header=FALSE)
colnames(eur) <- c("SNP_ID","ID","GENE_ID","PROBE_ID","CHR_SNP","CHR_GENE","SNPpos","TSSpos","distance","rvalue","pvalue")

afr <- fread(input=sprintf("zcat %s",afrgeu), sep='\t', header=FALSE)
colnames(afr) <- c("SNP_ID","ID","GENE_ID","PROBE_ID","CHR_SNP","CHR_GENE","SNPpos","TSSpos","distance","rvalue","pvalue")

merged <- merge(eur, afr, by="GENE_ID")
print(nrow(merged))

# color eQTL genes with significantly different effect sizes 
mean <- mean(abs(c(merged[,rvalue.x], merged[,rvalue.y])))
merged[, label := as.integer(abs(rvalue.x - rvalue.y) >= mean) ]

pdf("~/gtex-admix/plots/geu-eqtl-effect-sizes.pdf", width=12, height=6)
ggplot(merged, aes(x="rvalue.x",y="rvalue.y",colour=factor(label))) +
	geom_point() +
	labs(x="EUR eQTL r-value",y="AFR eQTL r-value") +
	theme_bw() +
	scale_colour_manual(values=c(1="red",0="black")) +
	theme(legend.position="none")
dev.off()





