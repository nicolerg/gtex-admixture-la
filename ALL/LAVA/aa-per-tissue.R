# get number of admixed individuals in each tissue
library(data.table)
library(ggplot2)

#pheno <- "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
# pheno <- fread(pheno, sep='\t', header=TRUE)
# # aa <- pheno[RACE==2|RACE==1,SUBJID] # AA or Asian
# aa <- pheno[RACE==2,SUBJID] # AA or Asian

admixed <- fread('/mnt/lab_data/montgomery/nicolerg/rfmix-admixed/combined/gtex-admixed0.9.txt', sep='\t', header=FALSE)
aa <- admixed[,V1]

flist <- dir(path = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices", 
	pattern = "normalized_expression.bed.gz", 
	full.names = TRUE)

t <- c()
n <- c()
frac <- c()
total <- c()
for (f in flist){
	if (!grepl("gz.tbi",f,fixed=TRUE)){
		x <- fread(cmd=sprintf("zcat %s",f),sep='\t',header=TRUE,nrows=3)
		tissue <- gsub(".*/","",f)
		tissue <- gsub(".v8.normalized_expression.bed.gz","",tissue,fixed=TRUE)
		tot <- ncol(x) - 4
		print(tissue)
		subjs <- colnames(x)
		subjs <- subjs[subjs %in% aa]
		naa <- length(subjs)
		print(naa)
		t <- c(t, tissue)
		n <- c(n, naa)
		total <- c(total, tot)
		print(naa/tot)
		frac <- c(frac, round(naa/tot,2))
	}
}

print(min(frac))
print(max(frac))

df <- data.frame(tissue=t,N_AA=n,frac=frac)
df$frac <- as.character(df$frac)
df$frac <- gsub("0.",".",df$frac,fixed=TRUE)

print(sum(total))
print(sum(n))

q()

write.table(df,'/mnt/lab_data/montgomery/nicolerg/local-eqtl/aa-per-tissue.tsv', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
#df <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/aa-per-tissue.tsv', sep='\t', header=TRUE)
#df$axiscol <- ifelse(df$N_AA < 57,'gray','black')

# tissuecols <- c(Adipose_Subcutaneous="#FF6600",
# 	Adipose_Visceral_Omentum="#FFAA00",
# 	Adrenal_Gland="#33DD33",
# 	Artery_Aorta="#FF5555",
# 	Artery_Coronary="#FFAA99",
# 	Artery_Tibial="#FF0000",
# 	Brain_Amygdala="#EEEE00",
# 	Brain_Anterior_cingulate_cortex_BA24="#EEEE00",
# 	Brain_Caudate_basal_ganglia="#EEEE00",
# 	Brain_Cerebellar_Hemisphere="#EEEE00",
# 	Brain_Cerebellum="#EEEE00",
# 	Brain_Cortex="#EEEE00",
# 	Brain_Frontal_Cortex_BA9="#EEEE00",
# 	Brain_Hippocampus="#EEEE00",
# 	Brain_Hypothalamus="#EEEE00",
# 	Brain_Nucleus_accumbens_basal_ganglia="#EEEE00",
# 	Brain_Putamen_basal_ganglia="#EEEE00",
# 	"Brain_Spinal_cord_cervical_c-1"="#EEEE00",
# 	Brain_Substantia_nigra="#EEEE00",
# 	Breast_Mammary_Tissue="#33CCCC",
# 	"Cells_EBV-transformed_lymphocytes"="#CC66FF",
# 	Cells_Cultured_fibroblasts="#AAEEFF",
# 	Colon_Sigmoid="#EEBB77",
# 	Colon_Transverse="#CC9955",
# 	Esophagus_Gastroesophageal_Junction="#8B7355",
# 	Esophagus_Mucosa="#552200",
# 	Esophagus_Muscularis="#BB9988",
# 	Heart_Atrial_Appendage="#9900FF",
# 	Heart_Left_Ventricle="#660099",
# 	Liver="#AABB66",
# 	Lung="#99FF00",
# 	Kidney_Cortex="yellow",
# 	Minor_Salivary_Gland="gray",
# 	Muscle_Skeletal="#AAAAFF",
# 	Nerve_Tibial="#FFD700",
# 	Ovary="#FFAAFF",
# 	Pancreas="#995522",
# 	Pituitary="#AAFF99",
# 	Prostate="#DDDDDD",
# 	Skin_Not_Sun_Exposed_Suprapubic="#0000FF",
# 	Skin_Sun_Exposed_Lower_leg="#7777FF",
# 	Small_Intestine_Terminal_Ileum="#555522",
# 	Spleen="#778855",
# 	Stomach="#FFDD99",
# 	Testis="#AAAAAA",
# 	Thyroid="#006600",
# 	Uterus="#FF66FF",
# 	Vagina="#FF5599",
# 	Whole_Blood="#FF00BB")

# labs <- c(Adipose_Subcutaneous="Subcutaneous adipose",
# 	Adipose_Visceral_Omentum="Visceral omentum",
# 	Adrenal_Gland="Adrenal gland",
# 	Artery_Aorta="Aorta",
# 	Artery_Coronary="Coronary artery",
# 	Artery_Tibial="Tibial artery",
# 	Brain_Amygdala="Amygdala",
# 	Brain_Anterior_cingulate_cortex_BA24="Anterior cingulate cortex",
# 	Brain_Caudate_basal_ganglia="Caudate nucleus (BG)",
# 	Brain_Cerebellar_Hemisphere="Cerebellar hemisphere",
# 	Brain_Cerebellum="Cerebellum",
# 	Brain_Cortex="Cortex",
# 	Brain_Frontal_Cortex_BA9="Frontal cortex",
# 	Brain_Hippocampus="Hippocampus",
# 	Brain_Hypothalamus="Hypothalamus",
# 	Brain_Nucleus_accumbens_basal_ganglia="Nucleus accumbens (BG)",
# 	Brain_Putamen_basal_ganglia="Putamen (BG)",
# 	"Brain_Spinal_cord_cervical_c-1"="Cervical spinal cord",
# 	Brain_Substantia_nigra="Substantia nigra",
# 	Breast_Mammary_Tissue="Breast mammary tissue",
# 	"Cells_EBV-transformed_lymphocytes"="EBV-transformed lymphocytes",
# 	Cells_Cultured_fibroblasts="Fibroblasts",
# 	Colon_Sigmoid="Sigmoid colon",
# 	Colon_Transverse="Transverse colon",
# 	Esophagus_Gastroesophageal_Junction="Gastroesophageal junction",
# 	Esophagus_Mucosa="Oesophagas mucosa",
# 	Esophagus_Muscularis="Oesophagas muscularis",
# 	Heart_Atrial_Appendage="Atrial appendage",
# 	Heart_Left_Ventricle="Left ventricle",
# 	Liver="Liver",
# 	Lung="Lung",
# 	Kidney_Cortex="Cortex",
# 	Minor_Salivary_Gland="Minor salivary gland",
# 	Muscle_Skeletal="Skeletal muscle",
# 	Nerve_Tibial="Tibial nerve",
# 	Ovary="Ovary",
# 	Pancreas="Pancreas",
# 	Pituitary="Pituitary",
# 	Prostate="Prostate",
# 	Skin_Not_Sun_Exposed_Suprapubic="NSE skin (suprapubic)",
# 	Skin_Sun_Exposed_Lower_leg="SE skin (lower leg)",
# 	Small_Intestine_Terminal_Ileum="Small intestine terminal ileum",
# 	Spleen="Spleen",
# 	Stomach="Stomach",
# 	Testis="Testis",
# 	Thyroid="Thyoid",
# 	Uterus="Uterus",
# 	Vagina="Vagina",
# 	Whole_Blood="Whole blood")

# pdf("~/gtex-admix/plots/figures/gtex-aa-pertissue.pdf", width=8.5, height=4)
# ggplot(df,aes(x=tissue,y=N_AA,fill=tissue)) +
# 	geom_bar(stat="identity") +
# 	geom_text(
# 		aes(label = frac, y = N_AA + 1),
# 		position = position_dodge(0.9),
# 		vjust = 0,
# 		size = 2) +
# 	scale_fill_manual(values=tissuecols,guide='none') +
# 	scale_x_discrete(labels=labs) +
# 	labs(x="GTEx v8 Tissue",y="N self-reported African Americans") +
# 	theme_bw() +
# 	theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size=8),
# 		axis.title.y=element_text(size=10),
# 		axis.title.x=element_blank()) 
# 	#geom_hline(yintercept=57,linetype="dashed")
# dev.off()
