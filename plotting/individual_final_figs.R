# final figures
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)
library(OneR)
library(ggrepel)
library(argparse)
library(plyr)
library(ggcorrplot)

master_data_dir <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/"
plot_dir <- "/users/nicolerg/gtex-admix/plots/new_peers"

parser <- ArgumentParser()
parser$add_argument("--admix_per_tissue", action="store_true")
parser$add_argument("--genotype_pc", action="store_true")
parser$add_argument("--local_block", action="store_true")
parser$add_argument("--local_block_trunc", action="store_true")
parser$add_argument("--egene_discovery", action="store_true")
parser$add_argument("--qq", action="store_true")
parser$add_argument("--pval_distn", action="store_true")
parser$add_argument("--asb_sub", action="store_true")
parser$add_argument("--venn", action="store_true")
parser$add_argument("--asb", action="store_true")
parser$add_argument("--tf_bind", action="store_true")
parser$add_argument("--cadd_plot", action="store_true")
parser$add_argument("--homer", action="store_true")
parser$add_argument("--local_var", action="store_true")
parser$add_argument("--tss_distance", action="store_true")
parser$add_argument("--coloc", action="store_true")
parser$add_argument("--coloc_fst", action="store_true")
parser$add_argument("--egene_fst", action="store_true")
parser$add_argument("--gtex_r2", action="store_true")
parser$add_argument("--fst_distn", action="store_true")
parser$add_argument("--max_fst_pop", action="store_true")
parser$add_argument("--coloc_locus_global", action="store_true")
parser$add_argument("--coloc_locus_local", action="store_true")
parser$add_argument("--gtex_esnps_r2_all", action="store_true")
parser$add_argument("--gtex_esnps_r2_high", action="store_true")
parser$add_argument("--local_var_with_fst", action="store_true")
parser$add_argument("--rfmix_pc_cor", action="store_true")
parser$add_argument("--cutoff", default=1e-06)

methodcols <- c(global='#81c0c7',LAVA='#CC7EAB')
methodlabs <- c(global='GlobalAC',LAVA='LocalAC')

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

shortlab <- c(Adipose_Subcutaneous='Subc.\nadipose',
	Artery_Tibial='Tibial\nartery',
	Lung='Lung',
	Muscle_Skeletal='Skeletal\nmuscle',
	Nerve_Tibial='Tibial\nnerve',
	Skin_Not_Sun_Exposed_Suprapubic='NSE skin')

short <- c(Adipose_Subcutaneous='Subc.\nadipose',
	Artery_Tibial='Artery',
	Lung='Lung',
	Muscle_Skeletal='Sk.\nmuscle',
	Nerve_Tibial='Nerve',
	Skin_Not_Sun_Exposed_Suprapubic='NSE\nskin')

tissuelabs <- c(Adipose_Subcutaneous="Subc. adipose",
	Artery_Tibial="Tibial artery",
	Lung="Lung",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Skin_Not_Sun_Exposed_Suprapubic="NSE skin")

mcnemar_exact <- function(a,b,c,d){
	
	x <- (b-c)^2/(b+c)
	pval <- pchisq(x, df=1, lower.tail=FALSE)

	return(pval)
}

tissuetextcol <- c(Adipose_Subcutaneous="#FF6600",
	Lung="#99FF00",
	Nerve_Tibial="#c6a700",
	Muscle_Skeletal="#AAAAFF",
	Skin_Not_Sun_Exposed_Suprapubic="#0000FF",
	Artery_Tibial="#FF0000")

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
	Minor_Salivary_Gland="#99BB88",
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

tissize <- c(Adipose_Subcutaneous=10,
	Adipose_Visceral_Omentum=8,
	Adrenal_Gland=8,
	Artery_Aorta=8,
	Artery_Coronary=8,
	Artery_Tibial=10,
	Brain_Amygdala=8,
	Brain_Anterior_cingulate_cortex_BA24=8,
	Brain_Caudate_basal_ganglia=8,
	Brain_Cerebellar_Hemisphere=8,
	Brain_Cerebellum=8,
	Brain_Cortex=8,
	Brain_Frontal_Cortex_BA9=8,
	Brain_Hippocampus=8,
	Brain_Hypothalamus=8,
	Brain_Nucleus_accumbens_basal_ganglia=8,
	Brain_Putamen_basal_ganglia=8,
	"Brain_Spinal_cord_cervical_c-1"=8,
	Brain_Substantia_nigra=8,
	Breast_Mammary_Tissue=8,
	"Cells_EBV-transformed_lymphocytes"=8,
	Cells_Cultured_fibroblasts=8,
	Colon_Sigmoid=8,
	Colon_Transverse=8,
	Esophagus_Gastroesophageal_Junction=8,
	Esophagus_Mucosa=8,
	Esophagus_Muscularis=8,
	Heart_Atrial_Appendage=8,
	Heart_Left_Ventricle=8,
	Liver=8,
	Lung=10,
	Kidney_Cortex=8,
	Minor_Salivary_Gland=8,
	Muscle_Skeletal=10,
	Nerve_Tibial=10,
	Ovary=8,
	Pancreas=8,
	Pituitary=8,
	Prostate=8,
	Skin_Not_Sun_Exposed_Suprapubic=10,
	Skin_Sun_Exposed_Lower_leg=8,
	Small_Intestine_Terminal_Ileum=8,
	Spleen=8,
	Stomach=8,
	Testis=8,
	Thyroid=8,
	Uterus=8,
	Vagina=8,
	Whole_Blood=8)

face <- c(Adipose_Subcutaneous='bold',
	Adipose_Visceral_Omentum='plain',
	Adrenal_Gland='plain',
	Artery_Aorta='plain',
	Artery_Coronary='plain',
	Artery_Tibial='bold',
	Brain_Amygdala='plain',
	Brain_Anterior_cingulate_cortex_BA24='plain',
	Brain_Caudate_basal_ganglia='plain',
	Brain_Cerebellar_Hemisphere='plain',
	Brain_Cerebellum='plain',
	Brain_Cortex='plain',
	Brain_Frontal_Cortex_BA9='plain',
	Brain_Hippocampus='plain',
	Brain_Hypothalamus='plain',
	Brain_Nucleus_accumbens_basal_ganglia='plain',
	Brain_Putamen_basal_ganglia='plain',
	"Brain_Spinal_cord_cervical_c-1"='plain',
	Brain_Substantia_nigra='plain',
	Breast_Mammary_Tissue='plain',
	"Cells_EBV-transformed_lymphocytes"='plain',
	Cells_Cultured_fibroblasts='plain',
	Colon_Sigmoid='plain',
	Colon_Transverse='plain',
	Esophagus_Gastroesophageal_Junction='plain',
	Esophagus_Mucosa='plain',
	Esophagus_Muscularis='plain',
	Heart_Atrial_Appendage='plain',
	Heart_Left_Ventricle='plain',
	Liver='plain',
	Lung='bold',
	Kidney_Cortex='plain',
	Minor_Salivary_Gland='plain',
	Muscle_Skeletal='bold',
	Nerve_Tibial='bold',
	Ovary='plain',
	Pancreas='plain',
	Pituitary='plain',
	Prostate='plain',
	Skin_Not_Sun_Exposed_Suprapubic='bold',
	Skin_Sun_Exposed_Lower_leg='plain',
	Small_Intestine_Terminal_Ileum='plain',
	Spleen='plain',
	Stomach='plain',
	Testis='plain',
	Thyroid='plain',
	Uterus='plain',
	Vagina='plain',
	Whole_Blood='plain')

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
	Kidney_Cortex="Cortex (kidney)",
	Minor_Salivary_Gland="Minor salivary gland",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Ovary="Ovary",
	Pancreas="Pancreas",
	Pituitary="Pituitary",
	Prostate="Prostate",
	Skin_Not_Sun_Exposed_Suprapubic="Skin (not sun-exposed)",
	Skin_Sun_Exposed_Lower_leg="Skin (sun-exposed)",
	Small_Intestine_Terminal_Ileum="Small intestine terminal ileum",
	Spleen="Spleen",
	Stomach="Stomach",
	Testis="Testis",
	Thyroid="Thyroid",
	Uterus="Uterus",
	Vagina="Vagina",
	Whole_Blood="Whole blood")

alltissues <- c('Adipose_Subcutaneous',
	'Adipose_Visceral_Omentum',
	'Adrenal_Gland',
	'Artery_Aorta',
	'Artery_Coronary',
	'Artery_Tibial',
	'Brain_Amygdala',
	'Brain_Anterior_cingulate_cortex_BA24',
	'Brain_Caudate_basal_ganglia',
	'Brain_Cerebellar_Hemisphere',
	'Brain_Cerebellum',
	'Brain_Cortex',
	'Brain_Frontal_Cortex_BA9',
	'Brain_Hippocampus',
	'Brain_Hypothalamus',
	'Brain_Nucleus_accumbens_basal_ganglia',
	'Brain_Putamen_basal_ganglia',
	"Brain_Spinal_cord_cervical_c-1",
	'Brain_Substantia_nigra',
	'Breast_Mammary_Tissue',
	"Cells_EBV-transformed_lymphocytes",
	'Cells_Cultured_fibroblasts',
	'Colon_Sigmoid',
	'Colon_Transverse',
	'Esophagus_Gastroesophageal_Junction',
	'Esophagus_Mucosa',
	'Esophagus_Muscularis',
	'Heart_Atrial_Appendage',
	'Heart_Left_Ventricle',
	"Liver",
	"Lung",
	'Kidney_Cortex',
	'Minor_Salivary_Gland',
	'Muscle_Skeletal',
	'Nerve_Tibial',
	"Ovary",
	"Pancreas",
	"Pituitary",
	"Prostate",
	'Skin_Not_Sun_Exposed_Suprapubic',
	'Skin_Sun_Exposed_Lower_leg',
	'Small_Intestine_Terminal_Ileum',
	"Spleen",
	"Stomach",
	"Testis",
	"Thyroid",
	"Uterus",
	"Vagina",
	'Whole_Blood')


admix_per_tissue <- function(w=4,h=6){
	
	aa <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/aa-per-tissue.tsv'
	aa <- fread(aa, sep='\t', header=TRUE)

	aa[,colour := ifelse(tissue %in% tissues, 'black', NA)]

	print(head(aa))
	
	g <- ggplot(aa,aes(x=tissue,y=N_AA,fill=tissue)) +
	geom_bar(stat="identity",aes(colour=colour)) +
	geom_text(
		aes(label = frac, y = N_AA + 2),
		#position = position_dodge(0.9),
		hjust = 0,
		vjust = 0.5,
		size = 2.5) +
	scale_fill_manual(values=tissuecols,guide='none') +
	scale_colour_manual(values=c(black='black','NA'=NA),guide='none') +
	scale_x_discrete(labels=labs,limits=rev(alltissues)) +
	labs(y="N admixed individuals (117AX)") +
	theme_bw() +
	theme(axis.text.y = element_text(hjust=1,vjust=0.5,size=rev(tissize),face=rev(face),colour='black'),
		axis.title.y=element_blank(),
		legend.margin=margin(t = 0, unit='mm'),
    	panel.grid.minor = element_blank(),
    	panel.grid.major.y = element_blank(),
    	panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_y_continuous(expand = c(0,0), limits=c(0,115)) +
    coord_flip()

	pdf("~/gtex-admix/plots/new_peers/admix_per_tissue.pdf",width=w,height=h)
	print(g)
	dev.off()

	return(g)
}




local_block <- function(w=6,h=4){

	# local ancestry tracks
	# read in global ancestry proportions
	indir <- "/mnt/lab_data/montgomery/nicolerg/rfmix-admixed/combined/"
	global_ai <- read.delim(paste0(indir,"gtex_ai_global_all.txt"), sep='\t', header=TRUE, stringsAsFactors=FALSE)

	# filter out subject with >0.9 homogeneity
	global_ai <- global_ai[ global_ai$AFR < 0.9 & global_ai$ASN < 0.9 & global_ai$EUR < 0.9 ,]
	admixed <- data.frame(ID=global_ai$ID)
	write.table(admixed,paste0(indir,"gtex-admixed0.9.txt"),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
	n <- nrow(admixed)
	# ID	ASN	EUR	AFR
	# GTEX-1117F	0.0064	0.1938	0.7997

	master <- data.frame(matrix(ncol = 6, nrow = 0))
	colnames(master) <- c("chr","start","stop","pop","ind","eur")

	inds <- global_ai$ID

	for (ind in inds){
		bedA <- paste0(indir,ind,"_A.bed")
		bedB <- paste0(indir,ind,"_B.bed")

		labs <- c("Chr","Start","Stop","Pop","Ind")
		df1 <- read.delim(bedA, sep='\t', header=FALSE, stringsAsFactors=FALSE)
		df1 <- df1[,1:4]
		df1$Ind <- paste(ind,"_A",sep='')
		colnames(df1) <- labs
		df2 <- read.delim(bedB, sep='\t', header=FALSE, stringsAsFactors=FALSE)
		df2 <- df2[,1:4]
		df2$Ind <- paste(ind,"_B",sep='')
		colnames(df2) <- labs

		df <- data.frame(rbind(df1, df2), stringsAsFactors=FALSE)
		df$Eur <- global_ai$EUR[match(ind, global_ai$ID)]

		# now add to master df
		master <- data.frame(rbind(master, df), stringsAsFactors=FALSE)
	}

	# standardize start positions and make a table of max stop values for each chromosome
	dict <- data.frame(c(1:22),rep(NA,22),rep(NA,22),rep(NA,22))
	colnames(dict) <- c("Chr","Max","Start","Stop")
	for (i in c(1:22)){
		rows <- which(master$Chr == i)
		min <- min(master$Start[rows])
		master$Start[rows] <- master$Start[rows] - min # make each chromosome start at 0
		master$Stop[rows] <- master$Stop[rows] - min # shift all segments over by that amount
		max <- max(master$Stop[rows])
		dict$Max[i] <- max
	}

	# figure out start and stop positions for each chromosome
	for (chr in c(1:22)){
		if (chr == 1){
			dict$Start[chr] <- 0
			dict$Stop[chr] <- dict$Max[chr]
		} else {
			dict$Start[chr] <- dict$Stop[chr-1] # make chromosomes start and stop on the same point
			dict$Stop[chr] <- dict$Start[chr] + dict$Max[chr]
		}
	}

	# add start value to every "start" and "stop" in master
	for (chr in c(1:22)){
		rows <- which(master$Chr == chr)
		master$Start[rows] <- master$Start[rows] + dict$Start[chr]
		master$Stop[rows] <- master$Stop[rows] + dict$Start[chr]
	}

	# sort master by descending eur values
	master <- master[order(master$Eur, decreasing=TRUE),] 
	# assign a y value to each haplotype (master$ind)
	master$Y <- NA
	indiv <- unique(master$Ind) #IND_A, #IND_B
	yval <- data.frame(indiv,c(1:length(indiv)),stringsAsFactors=FALSE)
	colnames(yval) <- c("ind","y")
	matchID <- function(x){
		yval$y[grep(x,yval$ind)]
	}
	master$Y <- unlist(lapply(master$Ind, matchID))
	master$Y <- as.numeric(master$Y)

	get_id <- function(x){
		return(unlist(unname(strsplit(x,'_')))[1])
	}
	ids <- unique(sapply(master$Ind, get_id))
	breaks <- seq(1.5,233.5,by=2)
	yint <- seq(2.5,232.5,by=2)

	g <- ggplot(master) + 
		geom_segment(aes(x=Start, xend=Stop, y=Y, yend=Y, colour=Pop), size=1, lineend="butt") +
		geom_hline(yintercept=yint, size=0.5) +
		geom_point(x=NA,y=NA,shape=22,aes(fill=Pop),colour='black',size=5) +
		theme_bw() +
		theme(panel.border = element_rect(size=1, colour='black'),
	          panel.grid.major = element_blank(),
	          panel.grid.minor = element_blank(), 
	          axis.ticks.x = element_blank(),
	          axis.ticks.y = element_blank(),
	          axis.text.y = element_blank(),
	          #axis.text.y = element_text(size = 4, colour='black'),
	          axis.text.x = element_text(size = 6, colour='black', angle=90, hjust=1, vjust=0.5),
	          plot.title = element_blank(),
	          #plot.title = element_text(hjust = 0.5, size = 14, colour='black'),
	          legend.text = element_text(colour='black',size=12),
	          legend.title = element_blank(),
	          legend.position = 'bottom',
	          axis.title.x = element_text(colour='black'),
	          axis.title.y = element_blank(),
	          legend.margin=margin(t = -3, unit='mm')) +
		scale_colour_manual(name="Local Ancestry", 
			values = c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
			guide = 'none') +
		scale_fill_manual(name="Local Ancestry", 
			values = c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
			breaks=c("ASN","EUR","AFR","UNK"),
			labels=c("Asian","European","African","Unknown")) +
		geom_vline(xintercept=dict$Stop[1:21], size=0.5) + 
		scale_x_discrete(limits=(dict$Start+dict$Stop)/2,
	        labels=c(1:22),name="Chromosome") +
		scale_y_continuous(breaks=breaks, labels=ids, expand = c(0,0), limits=c(0.5,234.5)) +
	  	guides(shape = guide_legend(override.aes = list(size = 4)))

	pdf("~/gtex-admix/plots/new_peers/block_haplo.pdf",width=w,height=h)
	print(g)
	dev.off()

	return(g)

}

rfmix_pc_cor <- function(w=6,h=2){
	
	pcs <- fread('~/gtex-admix/metadata/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',sep='\t',header=T)
	global_ai <- fread('~/gtex-admix/metadata/gtex_ai_global_all.txt',sep='\t',header=T)
	pcs[,SUBJID := sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
	m <- merge(pcs, global_ai, by.x='SUBJID', by.y='ID')
	m[,UNK := NULL]
	m <- data.frame(m)
	rownames(m) <- m$SUBJID
	m$SUBJID <- NULL
	m$FID <- NULL
	m$IID <- NULL
	corr <- cor(m)
	g <- ggcorrplot(corr[1:10, 21:23],lab=T, outline.col='black') + 
		theme(panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.text.x = element_text(hjust=1,colour='black'),
			axis.text.y = element_text(colour='black')) 

	pdf(paste0(plot_dir,"/rfmix_pc_cor.pdf"),width=w,height=h)
	print(g)
	dev.off()

	return(g)
}


egene_discovery <- function(h=3,w=5,cutoff=1e-6){

	# eGene discovery
	#load(paste0(master_data_dir,'egenes_master.RData'))
	load(paste0(master_data_dir,'egenes_master_fixed.RData'))
	egenes_master <- egenes_master_fixed 
	egenes_master <- egenes_master[pval_nominal < cutoff]

	counts <- data.table(TISSUE=tissues)
	counts[, N_diff_lead_snp := 0]
	counts[, N_same_lead_snp := 0] 
	counts[, N_global_only := 0]
	counts[, N_local_only := 0]
	counts[, total_global := 0]
	counts[, total_local := 0]

	for (t in tissues){
		
		sub <- egenes_master[tissue == t]
		sub[,pair := paste(gene_id, variant_id, sep=':')]
		global <- sub[method == 'global']
		local <- sub[method == 'LAVA']

		t_local <- length(unique(local[,gene_id]))
		t_global <- length(unique(global[,gene_id]))

		# called in both methods 
		g <- global[gene_id %in% local[,gene_id]]
		l <- local[gene_id %in% global[,gene_id]]

		repeats <- unique(g[D == 1, gene_id])
		#repeats <- unique(g[pair %in% l[,pair], gene_id])
		n_same_lead_snp <- length(repeats)
		n_diff_lead_snp <- length(unique(g[,gene_id])) - n_same_lead_snp

		global_only <- length(unique(global[!gene_id %in% local[,gene_id], gene_id]))
		local_only <- length(unique(local[!gene_id %in% global[,gene_id], gene_id]))

		counts[TISSUE == t, N_diff_lead_snp := n_diff_lead_snp]
		counts[TISSUE == t, N_same_lead_snp := n_same_lead_snp] 
		counts[TISSUE == t, N_global_only := global_only]
		counts[TISSUE == t, N_local_only := local_only]
		counts[TISSUE == t, total_global := t_global]
		counts[TISSUE == t, total_local := t_local]

	}

	egene_discovery <- ggplot(counts, aes(x=total_global, y=total_local, fill=factor(TISSUE))) +
		geom_point(colour='black', shape=21, size=2) +
		geom_abline(slope=1, intercept=0, linetype='dashed') +
		scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue', limits=tissues) +
		labs(x='N GlobalAC eGenes', y='N LocalAC eGenes') +
		theme_bw() +
		theme(#legend.text=element_text(size=10),
			#legend.title=element_text(size=12),
			legend.background = element_rect(color = "black", size = 0.5, linetype = "solid"),
			#legend.position = c(0.8,0.4),
			panel.grid.minor = element_blank(),
        	panel.border = element_rect(colour = "black", fill=NA, size=1)) 

	pdf(paste0(plot_dir,"/egene_discovery_",cutoff,".pdf"),width=w,height=h)
	print(egene_discovery)
	dev.off()

	return(egene_discovery)

}


qq <- function(w=5,h=4){

	if (!file.exists(paste0(master_data_dir,'qq.RData'))){

		load(paste0(master_data_dir,'egenes_master.RData'))
		load(paste0(master_data_dir,'allp_master.RData'))

		egenes_master[,pval_nominal := as.numeric(pval_nominal)]

		dtlist <- list()
		i <- 1
		for (t in tissues){
			sub <- allp_master[tissue == t]
			d <- as.data.table(qqplot(sub[,-log10(global_pval_nominal)], sub[,-log10(lava_pval_nominal)], plot.it=FALSE))
			colnames(d) <- c('global','lava')
			d[,tissue:=t]
			dtlist[[i]] <- d 
			i <- i+1 
		}

		qq <- data.table(rbindlist(dtlist))

		save(qq, file=paste0(master_data_dir,'qq.RData'))

	}else{
		load(paste0(master_data_dir,'qq.RData'))
	}

	q <- ggplot(qq, aes(x=global, y=lava, colour=factor(tissue))) +
		geom_line() +
		scale_colour_manual(values=tissuecols, labels=tissuelabs, name='GTEx v8 tissue', guide='none') +
		labs(x='GlobalAC (-log10 p-value)', y='LocalAC (-log10 p-value)') +
		geom_abline(slope=1, intercept=0, linetype='dashed')  +
		theme_bw() +
		theme(panel.grid.minor = element_blank(),
       		panel.border = element_rect(colour = "black", fill=NA, size=1))

	# # too much info to save as PDF 
	# jpeg(paste0(plot_dir,"/qq.jpg"),res=300,quality=100,width=w,height=h,units='in')
	# print(q)
	# dev.off()

	return(q)

}


pval_distn <- function(cutoff=1e-6,w=4,h=3){

	# # distribution of differences in p-values (unique eGenes at 1e-6 only)

	load(paste0(master_data_dir, '/uniq_egene_',cutoff,'.RData'))
	# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/allp_master-05.RData')
	load(paste0(master_data_dir, '/egenes_master.RData'))

	# split into unique local and global
	global <- uniq_egene[method=='global']
	local <- uniq_egene[method=='LAVA']

	global <- merge(global, egenes_master[method=='LAVA'], by=c('gene_id','tissue'), suffixes=c('.global','.local'))
	local <- merge(local, egenes_master[method=='global'], by=c('gene_id','tissue'), suffixes=c('.local','.global'))

	global <- global[,.(gene_id, tissue, pval_nominal.local, pval_nominal.global)]
	global[,method:='global']
	global <- unique(global, by=c('gene_id','tissue'))
	local <- local[,.(gene_id, tissue, pval_nominal.local, pval_nominal.global)]
	local[,method:='LAVA']
	local <- unique(local, by=c('gene_id','tissue'))

	global[,pval_diff := -log10(pval_nominal.global) + log10(pval_nominal.local)]
	local[,pval_diff := -log10(pval_nominal.local) + log10(pval_nominal.global)]

	uniq_egene_pval_diff <- data.table(rbind(global, local))

	g <- ggplot() +
		geom_histogram(data=uniq_egene_pval_diff[method=='LAVA'], aes(pval_diff, fill='LAVA'), alpha=1) +
		geom_histogram(data=uniq_egene_pval_diff[method=='global'], aes(pval_diff, fill='global'), alpha=0.7) +
		scale_fill_manual(values=methodcols, labels=methodlabs, breaks=c('LAVA','global')) +
		theme_bw() +
		theme(panel.grid.minor=element_blank(),
			panel.border = element_rect(colour = "black", fill=NA, size=1),
			legend.background = element_rect(color = "black", size = 0.5, linetype = "solid"),
			legend.position = c(0.7,0.80),
			legend.title = element_blank()) +
		labs(x='Difference in p-values (-log10)',y='N unique eGenes') 

	pdf(paste0(plot_dir,"/uniq_egene_pval_distn_",cutoff,".pdf"),width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

asb_sub <- function(cutoff=1e-5,w=3,h=3){

	# output of gtex_v8_rsid_map.py, allele_specific_binding.py
	asb_map <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/gtex_snp_to_asb_map.tsv'
	asb_map <- fread(asb_map, header=FALSE, sep='\t')
	colnames(asb_map) <- c('variant_id','asb_rsid')
	asb_map <- asb_map[!is.na(asb_map[,asb_rsid])]

	load(paste0(master_data_dir,'diff_lead_',cutoff,'.RData')) #diff_lead

	master_asb <- merge(diff_lead, asb_map, by='variant_id', all.x=TRUE)
	master_asb[,annotation := ifelse(is.na(asb_rsid),0,1)]

	dfsame <- list()
	i <- 1
	for (t in tissues){

		dt <- master_asb[tissue==t]

		all <- dt[,list(bin=max(annotation)), by = .(method, gene_id)]
		g <- all[method=='global',.(gene_id,bin)]
		l <- all[method=='LAVA',.(gene_id,bin)]
		all <- merge(g, l, by='gene_id', all=TRUE, suffixes=c('.global','.local'))

		total <- nrow(all)
		g_hit <- sum(all[,bin.global])
		l_hit <- sum(all[,bin.local])

		df2 <- data.frame(tissue=t,
			method=c('global','LAVA'),
			hit=c(g_hit, l_hit),
			count=total)

		dfsame[[i]] <- df2

		i <- i + 1

	}

	master2 <- data.table(rbindlist(dfsame))
	master2 <- master2[,.(method, hit)]
	m <- data.table(global=master2[method=='global',hit], LAVA=master2[method=='LAVA',hit])

	g2 <- ggplot(master2, aes(x=factor(method),y=hit,fill=factor(method)),colour='black') +
		geom_boxplot() +
		scale_fill_manual(values=methodcols, labels=methodlabs, name='Ancestry') +
		labs(y='N lead SNPs with ASB sites') +
		theme_bw() + 
		theme(axis.title.x=element_blank(),
			axis.text.x=element_text(colour='black',size=9),
			legend.position='none',
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			panel.border = element_rect(colour = "black", fill=NA, size=1),
			plot.margin = margin(1, 0.5, 1, 0.5, "cm")) +
		scale_x_discrete(labels=c(LAVA='Local',global='Global'),limits=c('global','LAVA')) +
		ylim(NA, 8.5) +
		geom_segment(x=1, xend=2, y=7.5, yend=7.5) +
		annotate("text",x=1.5, y=8.2, label='p = 0.01')
	
	pval <- signif(t.test(m[,global],m[,LAVA])$p.value, digits=3)
	writeLines('ASB sub')
	print(pval)

	pdf(paste0(plot_dir,"/asb-same-",cutoff,"-",pval,".pdf"),width=w,height=h)
	print(g2)
	dev.off()

	return(g2)

}

venn <- function(h=3.5,w=6.5,cutoff=1e-6){

	# discovery per tissue 

	load(paste0(master_data_dir,'egenes_master_fixed.RData'))
	egenes_master <- egenes_master_fixed
	egenes_master <- egenes_master[pval_nominal < cutoff]

	counts <- data.table(TISSUE=tissues)
	counts[, N_diff_lead_snp := 0]
	counts[, N_same_lead_snp := 0] 
	counts[, N_global_only := 0]
	counts[, N_local_only := 0]
	counts[, total_global := 0]
	counts[, total_local := 0]

	for (t in tissues){
		
		sub <- egenes_master[tissue == t]
		sub[,pair := paste(gene_id, variant_id, sep=':')]
		global <- sub[method == 'global']
		local <- sub[method == 'LAVA']

		t_local <- length(unique(local[,gene_id]))
		t_global <- length(unique(global[,gene_id]))

		# called in both methods 
		g <- global[gene_id %in% local[,gene_id]]
		l <- local[gene_id %in% global[,gene_id]]

		repeats <- unique(g[D == 1, gene_id])
		#repeats <- unique(g[pair %in% l[,pair], gene_id])
		n_same_lead_snp <- length(repeats)
		n_diff_lead_snp <- length(unique(g[,gene_id])) - n_same_lead_snp

		global_only <- length(unique(global[!gene_id %in% local[,gene_id], gene_id]))
		local_only <- length(unique(local[!gene_id %in% global[,gene_id], gene_id]))

		counts[TISSUE == t, N_diff_lead_snp := n_diff_lead_snp]
		counts[TISSUE == t, N_same_lead_snp := n_same_lead_snp] 
		counts[TISSUE == t, N_global_only := global_only]
		counts[TISSUE == t, N_local_only := local_only]
		counts[TISSUE == t, total_global := t_global]
		counts[TISSUE == t, total_local := t_local]

	}

	dt <- counts[,.(TISSUE, N_diff_lead_snp, N_same_lead_snp, N_global_only, N_local_only)]

	#dt <- fread('~/gtex-admix/metadata/eqtl_summ.txt', sep='\t', header=TRUE)

	dt[,local_y := N_global_only + N_diff_lead_snp + N_same_lead_snp + 0.5*N_local_only]
	dt[,global_y := 0.5*N_global_only]

	dt <- melt(dt, id=c('TISSUE','local_y','global_y'), variable.name='type', value.name='count')
	dt[,type := factor(type, levels = c('N_local_only', 'N_same_lead_snp', 'N_diff_lead_snp' ,'N_global_only'))]

	dt[type == 'N_global_only', global_label := count]
	dt[type == 'N_local_only', local_label := count]

	g <- ggplot(dt) +
	         geom_bar(aes(x=TISSUE, y=count, fill=type), stat='identity', colour='black') +
	         geom_text(aes(x=TISSUE, y=global_y, label=global_label)) +
	         geom_text(aes(x=TISSUE, y=local_y, label=local_label)) +
	         scale_fill_manual(values=c(N_diff_lead_snp='#8C7AA9',
	                                      N_same_lead_snp='gray',
	                                      N_global_only=unname(methodcols['global']),
	                                      N_local_only=unname(methodcols['LAVA'])),
	                             labels=c(N_diff_lead_snp='  Same eGene,\n  different eVariant',
	                                      N_same_lead_snp='  Same eGene,\n  same eVariant',
	                                      N_global_only='  eGene in GlobalAC only ',
	                                      N_local_only='  eGene in LocalAC only'),
	                             breaks=c('N_local_only', 'N_global_only', 'N_same_lead_snp', 'N_diff_lead_snp')) +
	         theme_bw() +
	         scale_x_discrete(labels=shortlab,limits=tissues) +
	         theme(axis.title.x=element_blank(),
	               axis.text.x=element_text(colour='black'),
	               #legend.text=element_text(size=10),
	               legend.title=element_blank(),
	               legend.position=c(0.5,0.85),
	               panel.grid.minor=element_blank(),
	               panel.border = element_rect(colour = "black", fill=NA, size=1),
	               legend.background = element_rect(color = "black", size = 0.5, linetype = "solid")) +
	         labs(y='N eGenes') +
	         ylim(c(0,max(dt[,local_y])+0.5*(max(dt[,local_y])))) +
	         guides(fill=guide_legend(ncol=2))

	pdf(paste0(plot_dir,"/venn_",cutoff,".pdf"),width=w,height=h)
	print(g)
	dev.off()

	return(g)

}


local_var <- function(h=3.5,w=8,cutoff=1e-04){

	# variance explained by local ancestry
	# refer to r2_all_snps.R 

	for (string in c('eur_afr_af_diff','eur_asn_af_diff','afr_asn_af_diff')){


		load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/master-admixed-all-merged.RData')
		
		# apply cutoff
		master <- master[is.na(pval_nominal.global) | (pval_nominal.global < cutoff & pval_nominal.local < cutoff)]

		# bin difference in allele frequency
		m <- master[(!is.na(gene_id)&!is.na(EUR_AF.global)) | is.na(gene_id)]
		m[,bin := bin(get(string), nbins=5, method='length')]

		bins <- as.character(unique(m[,bin]))

		# remove rand not in lowest bin
		m <- m[!(!is.na(local_rsq.rand) & !bin == bins[1])]

		# include only 3000 random 
		m[!is.na(local_rsq.rand), bin := paste0(bin,'_rand')]
		rand <- m[grepl("rand",bin)]
		m <- m[!grepl("rand",bin)]
		subset_rand <- rand[1:nrow(rand) %in% sample(nrow(rand), 3000, replace=FALSE)]
		m <- data.table(rbind(m, subset_rand))

		bins <- c(paste0(bins[1],'_rand'), bins)
		bin_n <- data.table(table(m[,bin]))
		colnames(bin_n) <- c('bin','label')

		m <- m[,.(bin, local_rsq.local, local_rsq.global, local_rsq.rand)]
		print(max(m[,local_rsq.local], na.rm=TRUE))

		new_bin <- paste0(bins[1])

		for (b in bin_n[,bin]){

			if(b==new_bin){next}
			bin_n[bin == b, pvalue := signif(t.test(m[bin == b, local_rsq.local], m[bin == b, local_rsq.global])$p.value, digits=3) ]

		}

		m <- melt(m, id='bin')
		m <- m[complete.cases(m)]

		if (cutoff==1e-04){

			if(string == 'eur_afr_af_diff'){ 

				labels <- c("(-0.000818,0.164]_rand"="(0,0.16]",
					"(-0.000818,0.164]"="(0,0.16]",
					"(0.164,0.327]"="(0.16,0.33]",
					"(0.327,0.491]"="(0.33,0.49]",
					"(0.491,0.654]"="(0.49,0.65]",
					"(0.654,0.819]"="(0.65,0.82]")

				limits <- c("(-0.000818,0.164]_rand",
					"(-0.000818,0.164]",
					"(0.164,0.327]",
					"(0.327,0.491]",
					"(0.491,0.654]",
					"(0.654,0.819]")

				x_lab <- 'Absolute difference b/w EUR and AFR AF in GlobalAC eVariant'

			} else if (string == 'eur_asn_af_diff'){

				labels <- c("(-0.000771,0.154]_rand"="(0,0.15]",
					"(-0.000771,0.154]"="(0,0.15]",
					"(0.154,0.308]"="(0.15,0.31]",
					"(0.308,0.462]"="(0.31,0.46]",
					"(0.462,0.616]"="(0.46,0.62]",
					"(0.616,0.771]"="(0.62,0.77]")

				limits <- c("(-0.000771,0.154]_rand",
					"(-0.000771,0.154]",
					"(0.154,0.308]",
					"(0.308,0.462]",
					"(0.462,0.616]",
					"(0.616,0.771]")

				x_lab <- 'Absolute difference b/w EUR and EAS AF in GlobalAC eVariant'

			} else if (string == 'afr_asn_af_diff'){

				limits <- c("(-0.000828,0.166]_rand",
					"(-0.000828,0.166]",
					"(0.166,0.331]",
					"(0.331,0.497]",
					"(0.497,0.662]",
					"(0.662,0.829]")

				labels <- c("(-0.000828,0.166]_rand"="(0,0.17]",
					"(-0.000828,0.166]"="(0,0.17]",
					"(0.166,0.331]"="(0.17,0.33]",
					"(0.331,0.497]"="(0.33,0.50]",
					"(0.497,0.662]"="(0.50,0.66]",
					"(0.662,0.829]"="(0.66,0.83]")

				x_lab <- 'Absolute difference b/w AFR and EAS AF in GlobalAC eVariant'

			}
		} else {
			limits <- bins
			labels <- bins
			x_lab=string
		}

		m[,variable := factor(variable, levels=c('local_rsq.global','local_rsq.local','local_rsq.rand'))]

		g <- ggplot() +
			geom_boxplot(data=m, aes(x=bin, y=sqrt(value), fill=variable), show.legend=FALSE, colour='black') +
			geom_point(data=m, size = -1, aes(x=bin, y=value, fill=variable), colour='black') + 
			scale_fill_manual(values=c(local_rsq.global=unname(methodcols['global']),local_rsq.local=unname(methodcols['LAVA']),local_rsq.rand='gray'),
				labels=c(local_rsq.global='GlobalAC eVariant ',local_rsq.local='LocalAC eVariant ',local_rsq.rand='Random variant '),breaks=c('local_rsq.rand','local_rsq.global','local_rsq.local')) +
			labs(x=x_lab,y='Correlation with LA') +
			theme_bw() +
			geom_text(data=bin_n, aes(x=bin, y=1, label=label)) +
			theme(axis.text.x=element_text(size=8),
				legend.text=element_text(size=12),
				legend.title=element_blank(),
	 			panel.grid.major=element_blank(),
	 			legend.background=element_blank(),
	 			#legend.background = element_rect(color = "black", size = 0.5, linetype = "solid"),
	 			panel.border = element_rect(colour = "black", fill=NA, size=1),
	 			legend.direction='horizontal',
	 			legend.position=c(0.5,0.92)) +
			guides(fill = guide_legend(override.aes = list(shape=22,size=6))) +
			scale_x_discrete(limits=limits,labels=labels) +
			scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.25,0.50,0.75,1.0))

		pdf(paste0(plot_dir,'/la_corr_stratified_by_delta_',string,'_',cutoff,'.pdf'),width=w,height=h)
		print(g)
		dev.off()

		if(string == 'eur_afr_af_diff'){
			g_eur_afr <- g
		}

	}

	return(g_eur_afr)

}

la_r2 <- function(cutoff=1e-04,w=4,h=4){

	if(!file.exists('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/m-r2_by_fst_bin.RData')){
		eur_afr_fst <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/high_local_fst/EUR_AFR.weir.fst'
		fst <- fread(eur_afr_fst, sep='\t', header=T)
		r2 <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/all.r2.txt',sep='\t',header=FALSE)
		colnames(r2) <- c('variant_id','r2')
		r2[,CHROM := sapply(variant_id, function(x) gsub('chr','',unname(unlist(strsplit(x,'_')))[1]) ) ]
		r2[,POS := sapply(variant_id, function(x) unname(unlist(strsplit(x,'_')))[2] ) ]
		m <- merge(fst, r2, by=c('CHROM','POS'))
		m[,fst := abs(as.numeric(WEIR_AND_COCKERHAM_FST))]
		m[,r2 := as.numeric(r2)]
		m <- m[,.(variant_id,fst,r2)]
		m <- m[complete.cases(m)]
		load(paste0(master_data_dir,'diff_lead_',cutoff,'.RData')) 
		m <- merge(m, diff_lead, by='variant_id',all.x=TRUE)

		g <- m[method=='global']
		l <- m[method=='LAVA']
		g[,method:=NULL]
		l[,method:=NULL]

		l <- l[,list(r2 = max(r2), fst = fst[which.max(r2)], variant_id = variant_id[which.max(r2)]), by=c('gene_id','tissue','pval_nominal')]
		g <- g[,list(r2 = max(r2), fst = fst[which.max(r2)], variant_id = variant_id[which.max(r2)]), by=c('gene_id','tissue','pval_nominal')]

		m <- merge(g, l, by=c('gene_id','tissue'), suffixes=c('_global','_local'))

		save(m, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/m-r2_by_fst_bin.RData')
	} else {
		load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/m-r2_by_fst_bin.RData')
	}

	print(nrow(m))

	m <- m[pval_nominal_local < cutoff & pval_nominal_global < cutoff]

	ld <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/plink-r2/local-global-pairs-out.txt',sep='\t',header=FALSE)
	ld <- unique(ld)
	colnames(ld) <- c('variant_id_global','variant_id_local','D')
	m_ld <- merge(m, ld, by=c('variant_id_local','variant_id_global'))

	pdf(paste0(plot_dir,'/diff-lead-ld-',cutoff,'.pdf'),width=w,height=h)
	g <- ggplot(m_ld, aes(D)) +
		geom_histogram() +
		theme_bw() +
		labs(x="LD (D')",y='Frequency') +
		theme(panel.grid.major=element_blank()) 
	print(g)
	dev.off()

	print(nrow(m_ld))
	print(nrow(m_ld[D < 1]))
	print(table(m_ld[D<1, tissue]))

	g <- ggplot(data=m_ld[D < 1],aes(x=sqrt(r2_global),y=sqrt(r2_local))) +
		geom_point(aes(colour=tissue), show.legend=FALSE, alpha=0.7) +
		theme_bw() +
		theme(panel.grid.major=element_blank()) +
		geom_abline(linetype='dashed',colour='red') +
		labs(x='GlobalAC correlation with LA',y='LocalAC correlation with LA')

	pdf(paste0(plot_dir,'/diff-lead-low-ld-',cutoff,'.pdf'),width=w,height=h)
	print(g)
	dev.off()

	print(head(m_ld))
	write.table(m_ld, '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/het_fraction/diff-lead-pairs-by-tissue-1e-04.txt',sep='\t',col.names=T, row.names=F, quote=F)

	print(t.test(m[,r2_global],m[,r2_local])$p.value)

	g <- ggplot() +
		geom_point(data=m, aes(x=fst_global, y=sqrt(r2_global), colour=tissue), show.legend=FALSE, alpha=0.7) +
		theme_bw() +
		theme(panel.grid.major=element_blank()) +
		geom_abline(linetype='dashed',colour='red') +
		labs(title='global')

	l <- ggplot() +
		geom_point(data=m, aes(x=fst_local, y=sqrt(r2_local), colour=tissue), show.legend=FALSE, alpha=0.7) +
		theme_bw() +
		theme(panel.grid.major=element_blank()) +
		geom_abline(linetype='dashed',colour='red') +
		labs(title='local')

	pdf(paste0(plot_dir,'/la-corr-global-v-local-',cutoff,'.pdf'),width=w,height=h)
	print(g)
	print(l)
	dev.off()


}
#la_r2()

local_var_with_fst <- function(){

	# variance explained by local ancestry
	# refer to r2_all_snps.R 



	bins <- as.character(unique(m[,bin]))
	bin_n <- data.table(table(m[,bin]))
	colnames(bin_n) <- c('bin','label')

	# breaks <- c("(-0.000798,0.16]_rand",
	# 	"(-0.000798,0.16]",
	# 	"(0.16,0.319]",
	# 	"(0.319,0.479]",
	# 	"(0.479,0.639]",
	# 	"(0.639,0.799]")

	# labels <- c("(-0.000798,0.16]_rand"='(0,0.16]',
	# 	"(-0.000798,0.16]"='(0,0.16]',
	# 	"(0.16,0.319]"='(0.16,0.32]',
	# 	"(0.319,0.479]"='(0.32,0.48]',
	# 	"(0.479,0.639]"='(0.48,0.64]',
	# 	"(0.639,0.799]"='(0.64,0.80]')

	m <- m[,.(bin, r2_local, r2_global, r2_rand)]

	#new_bin <- "(-0.000798,0.16]_rand"
	new_bin <- "(-0.000805,0.162]_rand"

	for (b in bin_n[,bin]){
		if(b==new_bin){next}
		bin_n[bin == b, pvalue := signif(t.test(m[bin == b, r2_local], m[bin == b, r2_global])$p.value, digits=3) ]

	}

	# print(mean(m[bin=="(-0.000798,0.16]",r2_global]))
	# print(mean(m[bin=="(-0.000798,0.16]",r2_local]))

	# print(mean(m[bin=="(0.639,0.799]",r2_global]))
	# print(mean(m[bin=="(0.639,0.799]",r2_local]))

	# print(mean(m[bin=="(0.639,0.799]",r2_global] - m[bin=="(0.639,0.799]",r2_local]))

	# print(t.test(m[bin=="(-0.000798,0.16]",r2_global], m[bin=="(-0.000798,0.16]_rand",r2_rand])$p.value)
	# print(t.test(m[bin=="(-0.000798,0.16]",r2_local], m[bin=="(-0.000798,0.16]_rand",r2_rand])$p.value)

	print(bin_n)

	m <- melt(m, id='bin')
	m <- m[complete.cases(m)]

	m[,variable := factor(variable, levels=c('r2_global','r2_local','r2_rand'))]

	g <- ggplot() +
		geom_boxplot(data=m, aes(x=bin, y=sqrt(value), fill=variable), show.legend=FALSE, colour='black') +
		geom_point(data=m, size = -1, aes(x=bin, y=value, fill=variable), colour='black') + 
		scale_fill_manual(values=c(r2_global=unname(methodcols['global']),r2_local=unname(methodcols['LAVA']),r2_rand='gray'),
			labels=c(r2_global='GlobalAC eVariant ',r2_local='LocalAC eVariant ',r2_rand='Random variant '),breaks=c('r2_rand','r2_global','r2_local')) +
		labs(x='Absolute EUR vs. AFR Fst of LocalAC eVariant',y='Correlation with LA') +
		theme_bw() +
		geom_text(data=bin_n, aes(x=bin, y=1, label=label)) +
		theme(axis.text.x=element_text(size=8),
			legend.text=element_text(size=12),
			legend.title=element_blank(),
 			panel.grid.major=element_blank(),
 			legend.background=element_blank(),
 			#legend.background = element_rect(color = "black", size = 0.5, linetype = "solid"),
 			panel.border = element_rect(colour = "black", fill=NA, size=1),
 			legend.direction='horizontal',
 			legend.position=c(0.5,0.92)) +
		guides(fill = guide_legend(override.aes = list(shape=22,size=6))) +
		# scale_x_discrete(limits=breaks,labels=labels,breaks=breaks) +
		scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.25,0.50,0.75,1.0))

	pdf(paste0(plot_dir,'/la_corr_stratified_by_delta_EUR-vs-AFR_fst_local_',cutoff,'.pdf'),width=8,height=3.5)
	print(g)
	dev.off()

	return(g)

}

# local_var <- function(h=6,w=8,cutoff=1e-04){

# 	load("/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/la_summary-all-annotations-merged.RData")
# 	print(head(la_summary))

# 	la <- la_summary[pval_nominal.local < cutoff & pval_nominal.global < cutoff]
# 	la <- la[lead_snp_D < 1]

# 	print(t.test(la[,rsq_genotype_local_ancestry.global], la[,rsq_genotype_local_ancestry.local]))

# 	la <- la[,.(tissue, gene_id, rsq_genotype_local_ancestry.global, rsq_genotype_local_ancestry.local, EUR_AFR_FST.global, EUR_AFR_FST.local)]
# 	print(head(la))

# 	la_g <- la[,.(tissue, gene_id, rsq_genotype_local_ancestry.global, EUR_AFR_FST.global)]
# 	la_g[,method:='global']
# 	la_l <- la[,.(tissue, gene_id, rsq_genotype_local_ancestry.local, EUR_AFR_FST.local)]
# 	la_l[,method:='LAVA']
# 	colnames(la_g) <- c('tissue','gene_id','rsq','fst','method')
# 	colnames(la_l) <- c('tissue','gene_id','rsq','fst','method')
# 	l <- rbind(la_g, la_l)

# 	# randomize row order
# 	l <- l[sample(1:nrow(l), nrow(l), replace=F)]

# 	g <- ggplot(l, aes(x=fst, y=rsq, colour=factor(method))) +
# 		geom_abline(linetype='dashed') +
# 		geom_point(alpha=0.5) +
# 		theme_bw() +
# 		scale_colour_manual(values=methodcols,labels=methodlabs) +
# 		theme(legend.title=element_blank(),
# 			panel.border = element_rect(colour = "black", fill=NA, size=1)) +
# 		labs(x='Absolute EUR vs. AFR Fst of eVariant',y='Variance in eVariant genotype explained by LA')

# 	pdf(paste0(plot_dir,"/fst_v_local_r2_",cutoff,".pdf"),width=w,height=h)
# 	print(g)
# 	dev.off()

# 	return(g)

# }

coloc <- function(w=11,h=6){

	# output of compile_coloc.R
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master_coloc-1e-04.RData')

	master_coloc[,label := NULL]

	gene_names <- fread('~/gtex-admix/metadata/gene_name_map.tsv', sep='\t', header=TRUE)

	master_coloc[, gene_stable_id := sapply(feature, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
	master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)

	master_coloc[,label_global:= ifelse(clpp_h4_local < 0.60 & clpp_h4_global > 0.75, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
	master_coloc[,label_local:= ifelse(clpp_h4_local > 0.75 & clpp_h4_global < 0.60, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

	master_coloc[,label := ifelse(((clpp_h4_local < 0.60 & clpp_h4_global > 0.75) | (clpp_h4_local > 0.75 & clpp_h4_global < 0.60)), gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

	g <- ggplot() +
	
		geom_rect(data=data.frame(xmin=0.75,xmax=1,ymin=0,ymax=0.6), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=0.5, fill=unname(methodcols['global'])) + # global is better
		geom_rect(data=data.frame(xmin=0,xmax=0.6,ymin=0.75,ymax=1), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=0.5, fill=unname(methodcols['LAVA'])) + # local is better

		geom_rect(data=data.frame(xmin=0.75,xmax=1,ymin=0,ymax=0.6), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), colour=unname(methodcols['global']),fill=NA) + # global is better
		geom_rect(data=data.frame(xmin=0,xmax=0.6,ymin=0.75,ymax=1), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), colour=unname(methodcols['LAVA']),fill=NA) + # local is better
		
		geom_label_repel(data=master_coloc[!is.na(label_local)],
			nudge_x=1.8 - master_coloc[!is.na(label_local),clpp_h4_global],
			direction="y",
			hjust=1,
			ylim=c(0.25,1),
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local,fill=tissue),
			seed=1,
			alpha=0.4) + 

		geom_label_repel(data=master_coloc[!is.na(label_local)],
			nudge_x=1.8 - master_coloc[!is.na(label_local),clpp_h4_global],
			direction="y",
			hjust=1,
			ylim=c(0.25,1),
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local),
			seed=1,
			fill=NA) + 

		geom_label_repel(data=master_coloc[!is.na(label_global)],
			nudge_x=1.8 - master_coloc[!is.na(label_global),clpp_h4_global],
			nudge_y=-0.5,
			direction="y",
			hjust=1,
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global,fill=tissue),
			seed=1,
			alpha=0.4) + 

		geom_label_repel(data=master_coloc[!is.na(label_global)],
			nudge_x=1.8 - master_coloc[!is.na(label_global),clpp_h4_global],
			nudge_y=-0.5,
			direction="y",
			hjust=1,
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global),
			seed=1,
			fill=NA) + 

		geom_point(data=master_coloc, aes(x=clpp_h4_global, y=clpp_h4_local, colour=tissue),alpha=0.7) +
		geom_point(data=master_coloc[!is.na(label_local) | !is.na(label_global)],alpha=1,aes(colour=tissue,clpp_h4_global, y=clpp_h4_local)) +

		scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue') +
		scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +

		labs(x='GlobalAC colocalization probability', y='LocalAC colocalization probability') +
		theme_bw() +
		geom_segment(data=data.frame(x=0,xend=1,y=0,yend=1),aes(x=x,y=y,xend=xend,yend=yend),colour='black',linetype='dashed') +
		geom_segment(data=data.frame(x=0,xend=1,y=0,yend=0),aes(x=x,y=y,xend=xend,yend=yend),size=0.3) +
		geom_segment(data=data.frame(x=0,xend=1,y=1,yend=1),aes(x=x,y=y,xend=xend,yend=yend),size=0.3) +
		geom_segment(data=data.frame(x=0,xend=0,y=1,yend=0),aes(x=x,y=y,xend=xend,yend=yend),size=0.3) +
		geom_segment(data=data.frame(x=1,xend=1,y=0,yend=1),aes(x=x,y=y,xend=xend,yend=yend),size=0.3) +
		scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1.8), expand = c(0,0)) +
		scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1), expand = c(0,0)) +
		guides(colour = guide_legend(override.aes=list(alpha=1,shape=15,size=4),ncol=6)) + 
		theme(panel.grid=element_blank(),
			legend.position='top',
			legend.justification='left',
			legend.title=element_blank(),
			legend.margin=margin(b = -3, unit='mm'),
			legend.direction='horizontal',
			axis.title.x=element_text(hjust=0.2),
			legend.text=element_text(size=16),
			panel.border=element_blank(),
			axis.title=element_text(size=16)) +
		annotate("text",x=0.8,y=0.15,hjust=0,label='Only\nGlobalAC\ncoloc. is\nsignificant',size=5) +
		annotate("text",x=0.05,y=0.9,hjust=0,label='Only LocalAC colocalization\nis significant',size=5) 
	
	pdf(paste0(plot_dir,'/COLOC-local-global-p4-scatter-bytissue-labelled.pdf'), width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

gtex_r2 <- function(w=6,h=3){
	# plot_gtex_la_r2.R, r2_gtex_lead_snps.R
	gtex_r2 <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/unique_lead_snp_r2.txt',sep='\t',header=FALSE)

	g <- ggplot(gtex_r2, aes(sqrt(V2))) +
		geom_histogram(bins=30) +
		theme_bw() + 
		labs(x='Correlation b/w GTEx\neVariant and LA',y='Frequency') +
		scale_x_continuous(expand=c(0,0),limits=c(0,0.8)) +
		scale_y_continuous(expand=c(0,0),limits=c(0,50000))

	pdf('~/gtex-admix/plots/new_peers/gtex-la-r2.pdf',width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

fst_distn <- function(w=4,h=4){
	# compile_coloc.R
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_all.RData')

	print(head(master))
	print(nrow(master))

	master[grepl('AFR',which_max_fst),superpop := 'Africa']
	master[grepl('EUR',which_max_fst),superpop := 'Europe']

	eur_max <- round(max(master[grepl('EUR',which_max_fst),max_fst]),digits=2)
	afr_max <- round(max(master[grepl('AFR',which_max_fst),max_fst]),digits=2)

	labels <- data.table(label=c(eur_max,afr_max),superpop=c('Europe','Africa'))
	print(table(master[,superpop]))
	print(head(master))

	g <- ggplot() +
		geom_histogram(data=master, aes(max_fst, fill=superpop), bins=40, colour='black') +
		geom_label(data=labels, aes(label=paste0('Max Fst: ',label)),x=Inf,y=Inf,hjust=1.1,vjust=1.2,label.size=NA,size=5) +
		theme_bw() + 
		labs(x='Max pairwise Fst value',y='Frequency') +
		scale_fill_manual(name='Superpopulation',values=c(Europe="#0000FF",Africa="#FF9900"),guide='none') +
		facet_wrap(~superpop, ncol=2) +
		scale_y_continuous(expand=c(0,0),limits=c(0,5.5e5)) +
		scale_x_continuous(expand=c(0,0),limits=c(0,0.2),breaks=c(0,0.05,0.1,0.15)) +
		theme(axis.title=element_text(size=16),
			strip.text = element_text(size=16))

	pdf('~/gtex-admix/plots/new_peers/fst_distn-merged-by-superpop.pdf', width=w, height=h)
	print(g)
	dev.off()

	return(g)
}

max_fst_pop <- function(w=5,h=4){
	# compile_coloc.R
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_all.RData')

	table <- data.table(table(master[,which_max_fst]))
	print(table)
	table[,N := N/nrow(master)*100]
	table <- table[order(N, decreasing=FALSE)]
	table[,colour := ifelse(grepl('EUR',V1), '#0000FF', '#FF9900')]
	table[,V1 := gsub('AFR_','',V1)]
	table[,V1 := gsub('EUR_','',V1)]
	order_of_labels <- table[,V1]

	g <- ggplot() +
		geom_bar(data=table, aes(x=V1, y=N, fill=colour), stat='identity', colour='black') +
		theme_bw() + 
		labs(x='Pairwise subpopulations',y='Frequency of max Fst value (%)') +
		scale_x_discrete(breaks=order_of_labels, limits=order_of_labels) +
		scale_fill_identity(guide='legend',name='Superpopulation',labels=c('#0000FF'='Europe','#FF9900'='Africa')) +
		theme(legend.position=c(0.6,0.13),
			legend.box.background = element_rect(colour = "black"),
			axis.text.y=element_text(colour='black',size=14),
			axis.title=element_text(size=16),
			legend.text=element_text(size=14),
			legend.title=element_text(size=16)) +
		coord_flip() +
		scale_y_continuous(expand = c(0,0),limits=c(0,max(table[,N]+0.5)))

	pdf('~/gtex-admix/plots/new_peers/fst_by_pair-all.pdf', width=w, height=h)
	print(g)
	dev.off()

	return(g)
}

coloc_locus <- function(method){
	# locus_plots.R
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master-locus_plots.RData')
	if (method == 'Global'){
		master <- master[global_p4 > local_p4]
		w = 9
		h = 8
		nc = 2
	} else if (method=='Local'){
		master <- master[local_p4 > global_p4]
		w = 17
		h = 12
		nc = 4
	} else {
		print('Method not recognized')
		return()
	}

	g <- ggplot(master) +
		geom_point(data=master[!is.na(variant_id)],colour='black',size=2,aes(x=pvalue_gwas,y=pvalue_eqtl)) +
		geom_point(aes(colour=method,x=pvalue_gwas,y=pvalue_eqtl)) +
		scale_colour_manual(values=c(Global=unname(methodcols['global']),Local=unname(methodcols['LAVA'])),labels=c(Global=unname(methodlabs['global']),Local=unname(methodlabs['LAVA']))) +
		theme_bw() +
		labs(x='GWAS p-value (-log10)',y='eQTL p-value (-log10)') +
		facet_wrap(~description,scales='free',ncol=nc) +
		theme(legend.title=element_blank())
	g <- g + geom_text(data = unique(master,by=c('description','p4_lab')), mapping = aes(x = -Inf, y = Inf, label = p4_lab),hjust = -0.1,vjust = 1.2)

	pdf(sprintf('~/gtex-admix/plots/new_peers/coloc-locus-compare-wrapped-%s.pdf',method),width=w,height=h)
	print(g)
	dev.off()

	return(g)

}

coloc_locus_global <- function(){
	g <- coloc_locus('Global')
	return(g)
}

coloc_locus_local <- function(){
	g <- coloc_locus('Local')
	return(g)
}

gtex_esnps_r2_all <- function(){
	
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2-all_signif_pairs.RData')
	gtex_esnp_r2[,pval_nominal := -log10(pval_nominal)]
	gtex_esnp_r2 <- gtex_esnp_r2[sample(1:nrow(gtex_esnp_r2))]
	
	g <- ggplot(gtex_esnp_r2, aes(colour=tissue, x=pval_nominal, y=sqrt(r2))) +
		geom_point(alpha=0.5) +
		theme_bw() +
		labs(x='GTEx v8 eQTL p-value (-log10)', y='Correlation b/w GTEx v8 eVariant and LA') +
		theme(legend.position='none') +
		scale_colour_manual(values=tissuecols) +
		scale_y_continuous(expand=c(0,0))

	return(g)

}

gtex_esnps_r2_high <- function(){
	
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2-all_signif_pairs.RData')
	gtex_esnp_r2[,pval_nominal := -log10(pval_nominal)]
	gtex_esnp_r2 <- gtex_esnp_r2[sample(1:nrow(gtex_esnp_r2))]
	x <- gtex_esnp_r2[r2>0.7]

	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_coloc-r2_all_esnps.RData')
	gtex_coloc <- gtex_coloc[p4>0.7]

	m <- merge(x, gtex_coloc, by=c('variant_id','gene_id','pval_nominal','r2','tissue'), all.x=TRUE)

	gene_names <- fread('~/gtex-admix/metadata/gene_name_map.tsv', sep='\t', header=TRUE)
	m[,gene_stable_id := sapply(gene_id, function(x) gsub("\\..*","",x))]
	m <- merge(m, gene_names, by='gene_stable_id')

	labs <- c(Skin_Sun_Exposed_Lower_leg='SE skin',
		Pancreas='Pancreas',
		Adipose_Visceral_Omentum='Visc. adipose',
		Whole_Blood='Blood',
		Artery_Tibial='Tibial artery')

	g <- ggplot(data=m, aes(x=pval_nominal, y=sqrt(r2))) +
		geom_point(data=m,alpha=0.7,aes(colour=tissue)) +
		theme_bw() +
		labs(x='GTEx eQTL p-value (-log10)', y='Correlation b/w GTEx eVariant and LA') +
		scale_colour_manual(values=tissuecols,breaks=unique(m[!is.na(trait),tissue]),limits=unique(m[!is.na(trait),tissue]),labels=labs) +
		geom_label_repel(data=m[!is.na(trait)],
			aes(label=gsub('_',' ',paste0(gene_name,': ',trait)),
				fill=tissue),
			nudge_x=37-m[!is.na(trait),pval_nominal],
			nudge_y=0.15,
			direction='y',
			hjust=0,
			seed=1,
			alpha=0.4,
			size=3,
			ylim=c(0.76,NA)) +
		geom_label_repel(data=m[!is.na(trait)],
			aes(label=gsub('_',' ',paste0(gene_name,': ',trait))),
			nudge_x=37-m[!is.na(trait),pval_nominal],
			nudge_y=0.15,
			direction='y',
			hjust=0,
			seed=1,
			fill=NA,
			alpha=1,
			size=3,
			ylim=c(0.76,NA)) +
		theme(legend.position='top',
			legend.title=element_blank(),
			legend.text=element_text(size=10),
			legend.margin=margin(b = -3, unit='mm'),) +
		guides(colour = guide_legend(override.aes=list(shape=15,size=4),ncol=5),fill=FALSE) +
		scale_fill_manual(values=tissuecols) +
		scale_x_continuous(expand=c(0,0), limits=c(0,130))
	return(g)
}


args <- parser$parse_args()

c <- as.numeric(args$cutoff)

if(args$gtex_esnps_r2_high){
	g <- gtex_esnps_r2_high()
}

if(args$gtex_esnps_r2_all){
	g <- gtex_esnps_r2_all()
}

if (args$coloc_locus_global){
	g <- coloc_locus_global()
}

if (args$coloc_locus_local){
	g <- coloc_locus_local()
}

if (args$gtex_r2){
	g <- gtex_r2()
}

if (args$fst_distn){
	g <- fst_distn()
}

if (args$max_fst_pop){
	g <- max_fst_pop()
}

if (args$admix_per_tissue){
	g <- admix_per_tissue()
}

if (args$genotype_pc){
	g <- genotype_pc()
}

if (args$local_block){
	g <- local_block()
}

if (args$local_block_trunc){
	g <- local_block_trunc()
}

if (args$egene_discovery){
	g <- egene_discovery(cutoff=c)
}

if (args$qq){
	g <- qq()
}

if (args$pval_distn){
	g <- pval_distn(cutoff=c)
}

if (args$asb_sub){
	g <- asb_sub(cutoff=c)
}

if (args$asb){
	g <- asb(cutoff=c)
}

if (args$venn){
	g <- venn(cutoff=c)
}

if (args$tf_bind){
	g <- tf_bind(cutoff=c)
}

if (args$cadd_plot){
	g <- cadd_plot(cutoff=c)
}

if (args$homer){
	g <- homer()
}

if (args$local_var){
	g <- local_var(cutoff=c)
}

if (args$coloc){
	g <- coloc()
}

if (args$tss_distance){
	g <- tss_distance(cutoff=c)
}

if (args$coloc_fst){
	g <- coloc_fst()
}

if (args$rfmix_pc_cor){
	g <- rfmix_pc_cor()
}

if (args$egene_fst){
	g <- egene_fst(cutoff=c)
}

if (args$local_var_with_fst){
	g <- local_var_with_fst()
}



# DEPRECATED ====================================================================

asb <- function(cutoff=1e-5,w=4,h=3){

	# output of gtex_v8_rsid_map.py, allele_specific_binding.py
	asb_map <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/gtex_mac10_snp_to_asb_map.tsv'
	asb_map <- fread(asb_map, header=FALSE, sep='\t')
	colnames(asb_map) <- c('variant_id','asb_rsid')
	asb_map <- asb_map[!is.na(asb_map[,asb_rsid])]

	load(paste0(master_data_dir,'diff_lead_',cutoff,'.RData')) 

	master_asb <- merge(diff_lead, asb_map, by='variant_id', all.x=TRUE)
	master_asb[,annotation := ifelse(is.na(asb_rsid),0,1)]

	dfsame <- list()
	counts_df <- list()
	summary_df <- list()
	all_df <- list()
	i <- 1
	for (t in tissues){

		dt <- master_asb[tissue==t]
		s <- data.frame(table(dt[,method],dt[,annotation]))
		s$tissue = t
		summary_df[[i]] <- s

		all <- dt[,list(bin=max(annotation)), by = .(method, gene_id)]
		g <- all[method=='global',.(gene_id,bin)]
		l <- all[method=='LAVA',.(gene_id,bin)]
		all <- merge(g, l, by='gene_id', all=TRUE, suffixes=c('.global','.local'))

		counts_df[[i]] <- all

		total <- nrow(all)
		g_hit <- sum(all[,bin.global])
		l_hit <- sum(all[,bin.local])

		df2 <- data.frame(tissue=t,
			method=c('global','LAVA'),
			hit=c(g_hit, l_hit),
			count=total)

		dfsame[[i]] <- df2

		# # get tissue-level p-values 
		# tab <- table(all[,bin.global], all[,bin.local])
		# print(t)
		# print(tab)
		# a <- tab[1,1]
		# b <- tab[1,2]
		# c <- tab[2,1]
		# d <- tab[2,2]
		# pval <- mcnemar_exact(a,b,c,d)
		# print(pval)

		i <- i + 1

		# print(sprintf('ASB at %s: %s: %s',cutoff, t, pval))

	}

	master2 <- data.table(rbindlist(dfsame))
	labels <- master2[,.(tissue, count)]
	labels <- unique(labels)

	label_pos <- 1.1*max(master2[,hit])

	#print(master2)

	if (cutoff==1e-06){
		
		lab = data.frame(tissue=c('Adipose_Subcutaneous',
									'Nerve_Tibial',
									'Artery_Tibial',
									'Muscle_Skeletal',
									'Lung',
									'Skin_Not_Sun_Exposed_Suprapubic'),
							method=c(rep('global',6),rep('LAVA',6)),
							label=c(0.58,0.95,0.97,0.90,0.27,1.0,
								0.90,0.53,0.61,0.04,0.78,0.24))
		lab = data.table(lab)
		# Adipose_Subcutaneous
		# global empirical p-value: 0.5802
		# local empirical p-value: 0.8986
		# Artery_Tibial
		# global empirical p-value: 0.9678
		# local empirical p-value: 0.611
		# Lung
		# global empirical p-value: 0.272
		# local empirical p-value: 0.7796
		# Muscle_Skeletal
		# global empirical p-value: 0.9046
		# local empirical p-value: 0.044
		# Nerve_Tibial
		# global empirical p-value: 0.9476
		# local empirical p-value: 0.5316
		# Skin_Not_Sun_Exposed_Suprapubic
		# global empirical p-value: 1.0
		# local empirical p-value: 0.2442
	} else if (cutoff==1e-05){

		lab = data.frame(tissue=c('Adipose_Subcutaneous',
									'Nerve_Tibial',
									'Artery_Tibial',
									'Muscle_Skeletal',
									'Lung',
									'Skin_Not_Sun_Exposed_Suprapubic'),
							method=c(rep('global',6),rep('LAVA',6)),
							label=c(0.39,0.95,0.74,0.33,0.43,0.82,
								0.55,0.53,0.69,0.01,0.34,0.20))
		lab = data.table(lab)
		# Adipose_Subcutaneous
		# global empirical p-value: 0.3946
		# local empirical p-value: 0.5492
		# Artery_Tibial
		# global empirical p-value: 0.7366
		# local empirical p-value: 0.6852
		# Lung
		# global empirical p-value: 0.429
		# local empirical p-value: 0.3354
		# Muscle_Skeletal
		# global empirical p-value: 0.3272
		# local empirical p-value: 0.012
		# Nerve_Tibial
		# global empirical p-value: 0.9476
		# local empirical p-value: 0.5334
		# Skin_Not_Sun_Exposed_Suprapubic
		# global empirical p-value: 0.8212
		# local empirical p-value: 0.1968
	}

	print(master2)
	master2 <- merge(master2, lab, by=c('tissue','method'))

	g2 <- ggplot() +
		geom_bar(data=master2, aes(x=factor(tissue),y=hit,fill=factor(method)), stat='identity', position='dodge',colour='black') +
		geom_text(data=labels, aes(label = count, x = tissue, y = label_pos), vjust = -0.5) +
		geom_text(data=master2, aes(x=factor(tissue),y=0.2,label=label,group=method), angle=90, hjust=0, position = position_dodge(0.9)) +
		scale_fill_manual(values=methodcols, labels=methodlabs, name='Ancestry') +
		labs(y='N lead SNPs with ASB sites') +
		theme_bw() + 
		theme(axis.title.x=element_blank(),
			axis.text.x=element_text(colour='black',size=9),
			legend.position='none',
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			panel.border = element_rect(colour = "black", fill=NA, size=1)) +
		scale_x_discrete(labels=short,limits=tissues) +
		ylim(0,1.2*max(master2[,hit])) 

	# # get pval across all tissues 
	counts <- data.table(rbindlist(counts_df))
	counts <- unique(counts)

	tab <- table(counts[,bin.global], counts[,bin.local])
	a <- tab[1,1]
	b <- tab[1,2]
	c <- tab[2,1]
	d <- tab[2,2]
	
	pval <- signif(mcnemar(a,b,c,d), digits=3)

	summary <- data.frame(rbindlist(summary_df))
	colnames(summary) <- c('method','hit','N','tissue')
	write.table(summary,paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/allele-specific-binding/summary_',cutoff,'.txt'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)

	pdf(paste0(plot_dir,"/asb-same-",cutoff,".pdf"),width=w,height=h)
	print(g2)
	dev.off()

	return(g2)
}
tf_bind <- function(cutoff=1e-5,h=3,w=4){

	# output of encode_tf_binding.py

	dfsame <- list()
	summary_df <- list()
	counts_df <- list()
	all_df <- list()
	i <- 1
	for (t in tissues){

		# define infiles 
		if (t == 'Muscle_Skeletal'){
			label <- 'gastrocnemius_medialis'
		} else if(t == 'Skin_Not_Sun_Exposed_Suprapubic'){
			label <- 'suprapubic_skin'
		} else if(t == 'Lung'){
			label <- 'lung'
		} else if(t == 'Nerve_Tibial'){
			label <- 'tibial_nerve'
		} else if(t == 'Artery_Tibial'){
			label <- 'tibial_artery'
		} else if(t == 'Adipose_Subcutaneous'){
			label <- 'subcutaneous_adipose_tissue'
		}

		infile <- paste0('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/encode-tf/',label,'/',label,'.tf.binding.tsv')
		dt <- fread(infile, sep='\t', header=FALSE)
		colnames(dt) <- c('gene_id','count','variant_id','pval_nominal','slope','method','annotation')
		dt <- dt[,.(gene_id, variant_id, method, annotation)]

		# same eGene, different lead SNP
		load(paste0(master_data_dir,'diff_lead_',cutoff,'.RData')) 
		diff_lead <- diff_lead[tissue == t]
		dt <- merge(dt, diff_lead, by=c('variant_id','gene_id','method'))

		s <- data.frame(table(dt[,method],dt[,annotation]))
		s$tissue = t
		summary_df[[i]] <- s

		all <- dt[,list(bin=max(annotation)), by = .(method, gene_id)]
		
		g <- all[method=='global',.(gene_id,bin)]
		l <- all[method=='LAVA',.(gene_id,bin)]
		all <- merge(g, l, by='gene_id', all=TRUE, suffixes=c('.global','.local'))

		all_df[[i]] <- all

		total <- nrow(all)
		g_hit <- sum(all[,bin.global])
		l_hit <- sum(all[,bin.local])

		df2 <- data.frame(tissue=t,
			method=c('global','LAVA'),
			hit=c(g_hit, l_hit),
			count=total)

		dfsame[[i]] <- df2

		# get tissue-level p-values 
		tab <- table(all[,bin.global], all[,bin.local])
		a <- tab[1,1]
		b <- tab[1,2]
		c <- tab[2,1]
		d <- tab[2,2]
		pval <- mcnemar_exact(a,b,c,d)

		print(sprintf('TF binding at %s: %s: %s',cutoff, t, pval))

		i <- i + 1

	}

	all <- rbindlist(all_df)
	all <- unique(all)
	print(head(all))
	tab <- table(all[,bin.global], all[,bin.local])
	print(tab)
	a <- tab[1,1]
	b <- tab[1,2]
	c <- tab[2,1]
	d <- tab[2,2]
	pval <- mcnemar_exact(a,b,c,d)
	print(pval)

	master2 <- data.table(rbindlist(dfsame))

	print(head(master2))

	label_pos <- 1.2*max(master2[,hit]) 

	g2 <- ggplot(master2, aes(x=factor(tissue),y=hit,fill=factor(method))) +
		geom_bar(stat='identity', position='dodge',colour='black') +
		geom_text(aes(label = count, x = tissue, y = 1.1*max(master2[,hit])), vjust = -0.5) +
		scale_fill_manual(values=methodcols, labels=methodlabs, name='Ancestry') +
		labs(y='N top SNP sets with TF binding sites') +
		theme_bw() +
		theme(axis.title.x=element_blank(),
 			axis.text.x=element_text(colour='black',size=9),
 			legend.position='none',
 			panel.grid.major=element_blank(),
 			panel.border = element_rect(colour = "black", fill=NA, size=1),
 			legend.background = element_rect(color = "black", size = 0.3, linetype = "solid")) +
		scale_x_discrete(labels=short,limits=tissues) +
		ylim(0,label_pos)

	# if (cutoff == 1e-5){
	# 	g2 <- g2 + geom_text(aes(x='Adipose_Subcutaneous',y=64,label='*'),colour='red',size=6)
	# }

	pdf(paste0(plot_dir,"/tf-binding-same-",cutoff,".pdf"),width=w,height=h)
	print(g2)
	dev.off()

	return(g2)
}
cadd_plot <- function(cutoff=1e-5,h=3,w=4){

	cadd <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/cadd/gtex.admixed.MAC10.all.snps_cadd.tsv'
	cadd <- fread(cadd, sep='\t', header=FALSE)
	colnames(cadd) <- c('variant_id','cadd_score')
	load(paste0(master_data_dir,'egenes_master.RData'))

	# same genes where lead SNP differs (diff_lead)

	load(paste0(master_data_dir,'diff_lead_',cutoff,'.RData'))
	cadd_master <- merge(diff_lead, cadd, by='variant_id', all.x=TRUE)

	c <- cadd_master[,.(gene_id, cadd_score, tissue, method)]
	c <- c[,list(max_score=max(cadd_score)), by = .(gene_id,tissue,method)]

	# get counts in each tissue 

	labels <- data.table(table(c[,tissue]))
	colnames(labels) <- c('tissue','count')
	labels[,count := count/2]

	labels[,y := max(c[,max_score],na.rm=TRUE)+0.2]

	pval <- signif(t.test(c[method=='global',max_score],c[method=='LAVA',max_score])$p.value,digits=3)

	writeLines('CADD scores')
	for (t in tissues){
		print(paste0(t,': ', t.test(c[tissue==t&method=='global',max_score],c[tissue==t&method=='LAVA',max_score])$p.value))
	}

	g2 <- ggplot() +
		geom_boxplot(data=c, aes(x=factor(tissue), y=max_score, fill=factor(method)),colour='black') +
		geom_text(data=labels, aes(x=factor(tissue), y=y, label=count),col='black') +
		labs(y='CADD score') + 
		theme_bw() +
		scale_fill_manual(values=methodcols, labels=methodlabs, name='Method') +
		theme(axis.title.x = element_blank(),
			axis.text.x=element_text(colour='black',size=9),
			legend.position='none',
			panel.grid.major=element_blank(),
			panel.border = element_rect(colour = "black", fill=NA, size=1),
 			legend.background = element_rect(color = "black", size = 0.3, linetype = "solid")) +
		scale_x_discrete(labels=short,limits=tissues) +
		ylim(c(min(c[,max_score],na.rm=TRUE),max(c[,max_score],na.rm=TRUE)+0.3))

	# if (cutoff == 1e-5){
	# 	g2 <- g2 + 
	# 		#geom_text(aes(x=c('Nerve_Tibial','Skin_Not_Sun_Exposed_Suprapubic'),y=3.3,label='*'),colour='red',size=6) +
	# 		geom_segment(aes(x=1.8,xend=2.2,y=3,yend=3)) +
	# 		geom_segment(aes(x=5.8,xend=6.2,y=3,yend=3)) +
	# 		annotate("text",
	# 			x = c(2,6),
	# 			y = c(3.3,3.3),
	# 			label = c("p = 0.002", "p = 0.01"))


	# }

	pdf(paste0(plot_dir,"/cadd-score-",cutoff,"-same-egenes-",pval,".pdf"),width=w,height=h)
	print(g2)
	dev.off()

	return(g2)
}
homer <- function(height=3,width=4){

	# refer to homer.R to see how this analysis was run 
	# tested: unique eQTLs in each method at pval < 0.01
	# regions: lead SNP + 15 bp on either side 
	# between global and local, one is the test region; the other is the background region

	# from looking at summary text file, FDR < 0.05:

	global <- c(0,0,0,0,1,0)
	local <- c(0,1,0,0,2,5)

	h <- data.table(tissue=tissues, global=global, LAVA=local)
	h <- melt(h, id='tissue')

	g <- ggplot(h, aes(x=factor(tissue), y=value, fill=factor(variable))) +
		geom_bar(stat='identity',position='dodge',colour='black') +
		scale_fill_manual(values=methodcols, labels=methodlabs, name='Ancestry') +
		labs(y='Number of enriched motifs (FDR < 0.05)') +
		theme_bw() +
		theme(axis.title.x=element_blank(),
				axis.text.x=element_text(colour='black',size=9),
				legend.position='none',
				panel.grid.major=element_blank(),
				panel.border = element_rect(colour = "black", fill=NA, size=1),
 				legend.background = element_rect(color = "black", size = 0.3, linetype = "solid")) +
		scale_x_discrete(labels=short,limits=tissues) +

	pdf("~/gtex-admix/plots/final_figures/homer-01-cutoff-all-diff-05FDR.pdf",width=width,height=height)
	print(g)
	dev.off()

	return(g)
}

coloc_fst <- function(){

	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master_coloc-1e-04.RData')

	master_coloc[,label := NULL]

	gene_names <- fread('~/gtex-admix/metadata/gene_name_map.tsv', sep='\t', header=TRUE)

	master_coloc[, gene_stable_id := sapply(feature, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
	master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)

	master_coloc[,label_global:= ifelse(clpp_h4_local < 0.60 & clpp_h4_global > 0.75, paste0(gene_name,':',gwas_trait), NA)]
	master_coloc[,label_local:= ifelse(clpp_h4_local > 0.75 & clpp_h4_global < 0.60, paste0(gene_name,':',gwas_trait), NA)]

	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-04.RData')
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_all.RData')

	global_better <- master_coloc[!is.na(label_global)]
	local_better <- master_coloc[!is.na(label_local)]

	global_better <- merge(global_better, diff_lead[method=='global'], by.y=c('tissue','gene_id'),by.x=c('tissue','feature'))
	global_better_fst <- merge(global_better, master, by='variant_id')

	# collapse down to single variant per gene with high Fst
	global_better_fst <- global_better_fst[,list(max_lead_fst=max(max_fst), max_fst_pop=which_max_fst[which.max(max_fst)]), by = .(tissue, feature)]

	local_better <- merge(local_better, diff_lead[method=='LAVA'], by.y=c('tissue','gene_id'),by.x=c('tissue','feature'))
	local_better_fst <- merge(local_better, master, by='variant_id')
	local_better_fst <- local_better_fst[,list(max_lead_fst=max(max_fst), max_fst_pop=which_max_fst[which.max(max_fst)]), by = .(tissue, feature)]

	g <- ggplot() +
		geom_histogram(data=master, aes(max_fst), bins=200) +
		geom_point(data=global_better_fst, aes(x=max_lead_fst,y=20000), fill=unname(methodcols['global']), shape=21, colour='black') +
		geom_point(data=local_better_fst, aes(x=max_lead_fst,y=50000), fill=unname(methodcols['LAVA']), shape=21, colour='black') +
		coord_cartesian(xlim=c(0,0.3)) +
		annotate('text',x=0.2,y=500000,label=paste0('Max Fst: ',max(master[,max_fst]))) +
		theme_bw() + 
		labs(x='Max pairwise Fst values',y='Count')

	pdf('~/gtex-admix/plots/new_peers/fst_distn-merged.pdf', width=6, height=3)
	print(g)
	dev.off()

	return(g)
}
egene_fst <- function(cutoff=c){

	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/diff_lead_1e-04.RData')
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_all.RData')
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/egenes_master.RData')
	egenes_master <- egenes_master[,.(gene_id,variant_id,pval_nominal,method,tissue)]
	g <- egenes_master[method=='global']
	l <- egenes_master[method=='LAVA']
	g[,method:=NULL]
	l[,method:=NULL]
	g[,pair := paste0(gene_id, ':', variant_id)]
	l[,pair := paste0(gene_id, ':', variant_id)]

	# remove egenes with overlapping lead SNPs
	tissues <- c('Adipose_Subcutaneous',
		'Nerve_Tibial',
		'Artery_Tibial',
		'Muscle_Skeletal',
		'Lung',
		'Skin_Not_Sun_Exposed_Suprapubic')

	diff_lead_global <- list()
	diff_lead_local <- list()
	i <- 1
	for (t in tissues){

		g_sub <- g[tissue == t]
		l_sub <- l[tissue == t]

		repeats <- l_sub[pair %in% g_sub[,pair], gene_id]
		l_sub <- l_sub[!(gene_id %in% repeats)]
		g_sub <- g_sub[!(gene_id %in% repeats)]

		diff_lead_global[[i]] <- g_sub
		diff_lead_local[[i]] <- l_sub

		i <- i + 1 

	}
	l <- data.table(rbindlist(diff_lead_local))
	g <- data.table(rbindlist(diff_lead_global))

	# I need to get Fst values for the missing SNPs (if they are in fact SNPs with multiple RSIDs for that position)

	nrow(egenes_master[!variant_id %in% master[,variant_id]])

	my_merge <- function(g, master){
		g <- merge(g, master, by='variant_id')
		g <- g[,list(max_fst=max(max_fst), variant_id=variant_id[which.max(max_fst)]), by = .(tissue, gene_id, pval_nominal)]
		return(g)
	}

	get_egenes <- function(master){
		g <- my_merge(g, master)
		l <- my_merge(l, master)

		egenes <- merge(g,l,by=c('gene_id','tissue'),suffixes=c('_local','_global'))
		egenes[,pval_diff := -log10(pval_nominal_local) + log10(pval_nominal_global)]

		egenes <- egenes[pval_nominal_local < cutoff | pval_nominal_global < cutoff]
		return(egenes)
	}

	egenes <- get_egenes(master)

	g <- ggplot(egenes) +
		geom_point(aes(x=max_fst_global, y=max_fst_local, colour=pval_diff), alpha=0.5) +
		scale_colour_gradient2(low=unname(methodcols['global']), high=unname(methodcols['LAVA']), mid='gray', name='-log10(Local p-val) -\n-log10(Global p-val)') +
		theme_bw() +
		labs(x='Max Global lead SNP Fst', y='Max Local lead SNP Fst') +
		geom_abline(linetype='dashed')
	
	pdf(paste0('~/gtex-admix/plots/new_peers/merged_fst_vs_lead_snp_pval_',cutoff,'.pdf'), width=5, height=3)
	print(g)
	dev.off()

	return(g)
}

tss_distance <- function(cutoff=1e-5,h=3,w=4){

	if(!file.exists(paste0(master_data_dir,'tss.RData'))){

		load(paste0(master_data_dir,'allp_master.RData'))
		load(paste0(master_data_dir,'egenes_master.RData'))

		x <- allp_master[, .(gene_id,variant_id,tissue,tss_distance) ]
		tss <- merge(egenes_master, x, by=c('gene_id','variant_id','tissue'))
		tss <- tss[!duplicated(tss)]
		tss[,pval_nominal := as.numeric(pval_nominal)]
		tss[,tss_distance := as.numeric(tss_distance)]
		save(tss,file=paste0(master_data_dir,'tss.RData'))

	} else {
		load(paste0(master_data_dir,'tss.RData'))
	}

	print(head(tss))

	# only same eGenes with different lead SNP 
	load(paste0(master_data_dir,'diff_lead_',cutoff,'.RData'))
	tss_master <- merge(tss, diff_lead, by=c('gene_id','variant_id','method','tissue'))
	# tss_master[tss_master == 1e7] <- NA

	tss_master <- tss_master[ , .SD[which.min(abs(tss_distance))], by = c('gene_id','method','tissue')]

	counts <- data.table(table(tss_master[method=='global',tissue]))
	colnames(counts) <- c('tissue','label')

	# counts <- as.data.table(table(dt[,variable]))
	# print(as.data.table(table(dt[,variable],dt[,tissue])))
	# lava_count <- counts[V1=='lava',N][1]
	# global_count <- counts[V1=='global',N][1]

	g1 <- ggplot() +
		geom_density(data=tss_master, alpha=0.2, aes(x=tss_distance, fill=method)) +
		scale_fill_manual(values=c(global='blue',LAVA='red'),name='Method',labels=c(global='Global',LAVA='Local')) + 
		labs(title='Same eGenes, different lead SNP',x='TSS distance',y='Density') +
		geom_vline(xintercept=0,linetype="dashed") +
		facet_wrap(~tissue,ncol=3) +
		geom_text(data=counts, aes(label=paste0('N = ',label), x=-750000, y=9e-7)) +
		theme_bw()
	
	pdf(paste0(plot_dir,"/tss_same_egenes_diff_lead_",cutoff,".pdf"),width=10,height=6)
	print(g1)
	dev.off()
}