library(data.table)
library(ggplot2)
outdir <- "~/gtex-admix/plots/figures/"

# FIGURE 1A. GTEx v8 genotype PCA =================================================================================================

pop <- fread('/users/nicolerg/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', header=TRUE)
pop <- pop[,list(SUBJID,RACE)]
pc <- fread('/users/nicolerg/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt', header=TRUE, sep='\t')
pc <- pc[,list(FID,PC1,PC2)]
split <- function(x){
	return(paste(unlist(unname(strsplit(x,'-')))[1:2], collapse='-'))
}
pc[,FID:=sapply(FID,split)]
m <- merge(pop, pc, by.x="SUBJID", by.y="FID")

admixed <- fread('/mnt/lab_data/montgomery/nicolerg/rfmix-admixed/combined/gtex-admixed0.9.txt',sep='\t',header=FALSE)
m[,shape:=ifelse(SUBJID %in% admixed[,V1],1,2)]

pdf(paste0(outdir,'GTEx_v8_genotype_pca.pdf'), height=4, width=4)
ggplot(data=NULL, aes(x=PC1,y=PC2)) +
	geom_point(data=m,aes(fill=factor(RACE),colour=factor(RACE),shape=factor(shape),alpha=factor(shape)),size=0.8) + 
	scale_fill_manual(values=c('1'="#ADFF33", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='gray'),
		labels=c('1'="Asian (N=11)", '2'="African American (N=103)", '3'="White (N=714)", '4'='American Indian (N=2)', '99'='Unknown (N=8)'),
		name="Self-Reported Race",guide='none') +
	scale_colour_manual(values=c('1'="#ADFF33", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='gray'),
		labels=c('1'="Asian (N=11)", '2'="African American (N=103)", '3'="White (N=714)", '4'='American Indian (N=2)', '99'='Unknown (N=8)'),
		name="Self-Reported Race",
		breaks=c('3','2','1','99','4')) +
	scale_shape_manual(name="Global Ancestry",labels=c('1'='<90% from single population','2'='>=90% from single population'),values=c('1'=21,'2'=16)) +
	scale_alpha_manual(values=c('1'=1,'2'=0.7), guide='none') +
	geom_point(data=m[shape==1],colour='black',shape=21,alpha=0.8) +
	theme_bw() +
	theme(legend.position=c(0.33,0.6),legend.background = element_rect(color = "black", size = 0.3, linetype = "solid")) +
	guides(colour = guide_legend(override.aes = list(size=2,alpha=1), order=1),
		shape = guide_legend(override.aes = list(size=c(2,1),alpha=1), order=0))
dev.off()



# FIGURE 1B. Admixed individuals per tissue =========================================================================================

aa <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/aa-per-tissue.tsv'
aa <- fread(aa, sep='\t', header=TRUE)

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

labs <- c(Adipose_Subcutaneous="Subcutaneous adipose",
	Adipose_Visceral_Omentum="Visceral omentum",
	Adrenal_Gland="Adrenal gland",
	Artery_Aorta="Aorta",
	Artery_Coronary="Coronary artery",
	Artery_Tibial="Tibial artery",
	Brain_Amygdala="Amygdala",
	Brain_Anterior_cingulate_cortex_BA24="Anterior cingulate cortex",
	Brain_Caudate_basal_ganglia="Caudate nucleus (BG)",
	Brain_Cerebellar_Hemisphere="Cerebellar hemisphere",
	Brain_Cerebellum="Cerebellum",
	Brain_Cortex="Cortex",
	Brain_Frontal_Cortex_BA9="Frontal cortex",
	Brain_Hippocampus="Hippocampus",
	Brain_Hypothalamus="Hypothalamus",
	Brain_Nucleus_accumbens_basal_ganglia="Nucleus accumbens (BG)",
	Brain_Putamen_basal_ganglia="Putamen (BG)",
	"Brain_Spinal_cord_cervical_c-1"="Cervical spinal cord",
	Brain_Substantia_nigra="Substantia nigra",
	Breast_Mammary_Tissue="Breast mammary tissue",
	"Cells_EBV-transformed_lymphocytes"="EBV-transformed lymphocytes",
	Cells_Cultured_fibroblasts="Fibroblasts",
	Colon_Sigmoid="Sigmoid colon",
	Colon_Transverse="Transverse colon",
	Esophagus_Gastroesophageal_Junction="Gastroesophageal junction",
	Esophagus_Mucosa="Oesophagas mucosa",
	Esophagus_Muscularis="Oesophagas muscularis",
	Heart_Atrial_Appendage="Atrial appendage",
	Heart_Left_Ventricle="Left ventricle",
	Liver="Liver",
	Lung="Lung",
	Kidney_Cortex="Cortex",
	Minor_Salivary_Gland="Minor salivary gland",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Ovary="Ovary",
	Pancreas="Pancreas",
	Pituitary="Pituitary",
	Prostate="Prostate",
	Skin_Not_Sun_Exposed_Suprapubic="NSE skin (suprapubic)",
	Skin_Sun_Exposed_Lower_leg="SE skin (lower leg)",
	Small_Intestine_Terminal_Ileum="Small intestine terminal ileum",
	Spleen="Spleen",
	Stomach="Stomach",
	Testis="Testis",
	Thyroid="Thyoid",
	Uterus="Uterus",
	Vagina="Vagina",
	Whole_Blood="Whole blood")

pdf("~/gtex-admix/plots/figures/gtex-admix-pertissue.pdf", width=8.5, height=4)
ggplot(aa,aes(x=tissue,y=N_AA,fill=tissue)) +
	geom_bar(stat="identity") +
	geom_text(
		aes(label = frac, y = N_AA + 1),
		position = position_dodge(0.9),
		vjust = 0,
		size = 2) +
	scale_fill_manual(values=tissuecols,guide='none') +
	scale_x_discrete(labels=labs) +
	labs(x="GTEx v8 Tissue",y="N admixed individuals") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size=8,colour='black'),
		axis.title.y=element_text(size=14),
		axis.title.x=element_blank()) 
	#geom_hline(yintercept=57,linetype="dashed")
dev.off()



# FIGURE 1C. Discovery for different eQTL calling methods ========================================================================

for (tissue in c("Lung","Thyroid")){
	globaleqtl <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/',tissue,'/gtex.admix.global.eqtl.txt.gz')
	localeqtl <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/',tissue,'/gtex.admix.lava.eqtl.txt.gz')

	ge <- fread(input=sprintf("zcat %s",globaleqtl),sep='\t',header=TRUE)
	le <- fread(input=sprintf("zcat %s",localeqtl),sep='\t',header=TRUE)

	n_lava_eqtl <- nrow(le)
	n_global_eqtl <- nrow(ge)
	m <- merge(ge,le,by="gene_id")
	n_shared_egenes <- nrow(m)
	n_shared_esnp <- nrow(m[variant_id.x == variant_id.y])

	df <- data.frame(x=c("LAVA","global","eGene","eSNP"),
		y=c(n_lava_eqtl,n_global_eqtl,n_shared_egenes,n_shared_esnp))
	df$tissue <- tissue

	assign(paste0(tissue,'_df'),df)
}

df <- data.frame(rbind(Lung_df, Throid_df))

colours <- c(Lung="#99FF00",Thyroid="#006600")

pdf(paste0(outdir,"eqtl-discoveries.pdf"),height=4,width=6)
ggplot(df, aes(x=x,y=y,fill=factor(tissue))) +
	geom_bar(stat="identity",colour='black') +
	labs(y="N eGenes") +
	theme_bw() +
	theme(axis.text.x=element_text(colour="black",size=12),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=12)) +
	scale_x_discrete(limits=c("LAVA","global","eGene","eSNP"),
		labels=c(LAVA="LAVA",global="global",eGene="shared\neGenes",eSNP="same\nlead SNP")) +
	scale_fill_manual(name='GTEx Tissue',values=colours)
dev.off()



# FIGURE 2A. EUR into AFR and AFR into EUR discoveries ========================================================================

m <- fread("/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/merged-eqtls.tsv", sep='\t', header=TRUE)

pdf(paste0(outdir,"geu-eqtl-by-pop.pdf"),width=4,height=4)
ggplot(data=NULL) +
	geom_point(data=m,aes(x=slope.EUR,y=slope.AFR,colour=factor(pop)),alpha=0.7,size=0.5) +
	labs(x="EUR effect size",y="AFR effect size") +
	scale_colour_manual(values=c(EUR2AFR="#0000FF",AFR2EUR="#FF9900"),
		name="eQTL population",
		labels=c(EUR2AFR='EUR into AFR', AFR2EUR='AFR into EUR'),
		guide='legend') +
	geom_point(data=m[sigpoint==1], aes(x=slope.EUR,y=slope.AFR,fill=factor(pop),shape=factor(sigpoint)), size=1, alpha=1, colour="black") +
	scale_shape_manual(values=21,labels='Top 100 most\ndifferent') +
	scale_fill_manual(name="Different effect sizes",
		values=c(EUR2AFR="#0000FF",AFR2EUR="#FF9900"),
		guide='none') +
	guides(colour = guide_legend(override.aes = list(size=2),order=1,keywidth=unit(0.9,"cm")),
		shape = guide_legend(override.aes = list(size=2),order=0,label.vjust=1.5,keywidth=unit(1,"cm"))) +
	# ylim(lims) +
	# xlim(lims) +
	theme_bw() +
	theme(legend.position=c(0.79,0.23),
		legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"),
		legend.title=element_blank())
dev.off()


