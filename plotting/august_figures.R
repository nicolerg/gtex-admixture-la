library(data.table)
library(ggplot2)
library(reshape)
library(OneR)
library(ggrepel)
library(argparse)
library(plyr)
library(ggcorrplot)
library(gridExtra)
library(gtable)
library(grid)

master_data_dir = "/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/merged/"
plot_dir = "/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/plots/"
supp_dir = "/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/plots/supplement/"
meta_dir = "/oak/stanford/groups/smontgom/nicolerg/LAVA/REVISIONS/metadata/"

parser <- ArgumentParser()
# individual figures 
parser$add_argument("--genotype_pc", action="store_true")
parser$add_argument("--admix_per_tissue_short", action="store_true")
parser$add_argument("--rfmix_pc_cor", action="store_true")
parser$add_argument("--local_block_trunc", action="store_true")
parser$add_argument("--qq", action="store_true")
parser$add_argument("--egene_discovery", action="store_true")
parser$add_argument("--venn", action="store_true")
parser$add_argument("--pval_distn", action="store_true")
parser$add_argument("--gtex_r2", action="store_true")
parser$add_argument("--gtex_esnps_r2_high", action="store_true")
parser$add_argument("--coloc", action="store_true")
parser$add_argument("--fst_distn", action="store_true")
parser$add_argument("--max_fst_pop", action="store_true")
parser$add_argument("--regression_res", action="store_true")
parser$add_argument("--ve_expr_ancestry", action="store_true")
# supplement
parser$add_argument("--admix_per_tissue", action="store_true")
parser$add_argument("--locus_plots", action="store_true")
parser$add_argument("--tg_gtex_pca", action="store_true")
# main figures
parser$add_argument("--figure1", action="store_true")
parser$add_argument("--figure2", action="store_true")
parser$add_argument("--figure3", action="store_true")
parser$add_argument("--figure4", action="store_true")
parser$add_argument("--cutoff", default=1e-06)

methodcols <- c(global='#81c0c7',LAVA='#CC7EAB')
methodlabs <- c(global='GlobalAA',LAVA='LocalAA')

tissues <- c('Nerve_Tibial',
	'Artery_Tibial',
	'Adipose_Subcutaneous',
	'Whole_Blood',
	'Muscle_Skeletal',
	'Skin_Not_Sun_Exposed_Suprapubic',
	'Lung'
	)

shortlab <- c(Adipose_Subcutaneous='Subc.\nadipose',
	Artery_Tibial='Tibial\nartery',
	Lung='Lung',
	Muscle_Skeletal='Skeletal\nmuscle',
	Nerve_Tibial='Tibial\nnerve',
	Skin_Not_Sun_Exposed_Suprapubic='NSE\nskin',
	Whole_Blood='Whole\nblood')

tissuelabs <- c(Adipose_Subcutaneous="Subc. adipose",
	Artery_Tibial="Tibial artery",
	Lung="Lung",
	Muscle_Skeletal="Skeletal muscle",
	Nerve_Tibial="Tibial nerve",
	Skin_Not_Sun_Exposed_Suprapubic="NSE skin",
	Whole_Blood='Whole blood')

tissuetextcol <- c(Adipose_Subcutaneous="#FF6600",
	Lung="#99FF00",
	Nerve_Tibial="#c6a700",
	Muscle_Skeletal="#AAAAFF",
	Skin_Not_Sun_Exposed_Suprapubic="#0000FF",
	Artery_Tibial="#FF0000",
	Whole_Blood="#FF00BB")

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
	Lung="green4", # changed from #9aff00, which is lime green
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
	Whole_Blood=10)

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
	Whole_Blood='bold')

all_labs <- c(Adipose_Subcutaneous="Subc. adipose",
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
	Skin_Not_Sun_Exposed_Suprapubic="NSE skin",
	Skin_Sun_Exposed_Lower_leg="SE skin",
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

genotype_pc <- function(w=4,h=4){

	# genotype PC plot with selected individuals shown 

	pop <- fread(sprintf('%s/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt',meta_dir), sep='\t', header=TRUE)
	pop <- pop[,list(SUBJID,RACE)]
	pc <- fread(sprintf('%s/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',meta_dir), header=TRUE, sep='\t')
	pc <- pc[,list(FID,PC1,PC2)]
	split <- function(x){
		return(paste(unlist(unname(strsplit(x,'-')))[1:2], collapse='-'))
	}
	pc[,FID:=sapply(FID,split)]
	m <- merge(pop, pc, by.x="SUBJID", by.y="FID")

	admixed <- fread(sprintf('%s/gtex-admixed0.9.txt',meta_dir),sep='\t',header=FALSE)
	m[,shape:=ifelse(SUBJID %in% admixed[,V1],1,2)]

	m <- data.table(rbind(m, data.table(SUBJID='filler',
		RACE=10,
		PC1=NA,
		PC2=NA,
		shape=1)))
	print(m)

	print(table(m[,RACE]))

	g <- ggplot(data=NULL, aes(x=PC1,y=PC2)) +
		geom_point(data=m,aes(fill=factor(RACE),colour=factor(RACE),shape=factor(shape),alpha=factor(shape),size=factor(shape))) + 
		scale_fill_manual(values=c('1'="green3", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='gray', '10'=NA),guide='none') +
		scale_colour_manual(values=c('1'="green3", '2'="#FF9900", '3'="#0000FF", '4'='#FF3D3D', '99'='gray','10'='black'),
			labels=c('1'="Asian American (N=12)", '2'="African American (N=103)", '3'="European American (N=714)", '4'='Native American (N=2)', '99'='Unknown (N=8)'),
			breaks=c('3','2','1','99','4')) +
		scale_size_manual(values=c('1'=2,'2'=1), guide='none') +
		scale_shape_manual(values=c('1'=21,'2'=16), breaks='2', limits='2', labels='117 admixed (117AX)') +
		scale_alpha_manual(values=c('1'=0.8,'2'=0.6), guide='none') +
		geom_point(data=m[shape==1],colour='black',shape=21,alpha=1,size=2,aes(fill=factor(RACE))) +
		geom_point(data=m[shape==1],colour='black',shape=21,alpha=1,size=2,fill=NA) +
		theme_classic() +
		theme(legend.justification=c(0,0), 
			legend.position=c(0.01,0.51),
			legend.text = element_text(size=10),
			legend.title = element_blank(),
        	legend.spacing = unit(3,unit='mm')) +
		guides(colour = guide_legend(override.aes = list(size=2,
									alpha=1,
									shape=21,
									colour=NA,
									fill=c("#0000FF","#FF9900","green3",'gray','#FF3D3D')), keywidth=5, keyheight=5, default.unit="mm"),
				shape = guide_legend(override.aes = list(size=2,
					alpha=1,
					shape=21,
					colour='black',fill='white'), keywidth=5, keyheight=5, default.unit="mm")) +
		labs(x='gPC1',y='gPC2')

	pdf(sprintf("%s/fig1a_geno_pc.pdf",plot_dir),width=w,height=h)
	print(g)
	dev.off()

	return(g)

}

admix_per_tissue_short <- function(w=4,h=4){

	aa <- sprintf('%s/aa-per-tissue.tsv',meta_dir)
	aa <- fread(aa, sep='\t', header=TRUE)

	print(head(aa))
	aa <- aa[order(N_AA,decreasing=T)]

	aa[,colour := ifelse(tissue %in% tissues, 'black', NA)]
	aa <- aa[N_AA > 30]
	limits <- rev(aa[,tissue])

	# change tissue colors
	tcols <- tissuecols
	tcols[!names(tcols) %in% tissues] <- 'gray'

	# print(nrow(aa))
	# sizes = unname(unlist(tissize[match(limits, names(tissize))]))
	# print(sizes)
	
	g <- ggplot(aa,aes(x=tissue,y=N_AA,fill=tissue)) +
		geom_bar(stat="identity",aes(colour=colour)) +
		scale_fill_manual(values=tcols,guide='none') +
		scale_colour_manual(values=c(black='black','NA'=NA),guide='none') +
		scale_x_discrete(labels=all_labs,limits=limits) +
		labs(y="N admixed individuals (117AX)") +
		theme_classic() +
		theme(axis.text.y = element_text(hjust=1,vjust=0.5,
			size=tissize[match(limits, names(tissize))],
			face=face[match(limits, names(face))],
			colour='black'),
			axis.title.y=element_blank()) +
	 	scale_y_continuous(expand = c(0,0), limits=c(0,115)) +
	    coord_flip()

    #p = g + theme(axis.text.y = element_text(size=sizes))

	pdf(sprintf("%s/fig1c_admix_per_tissue_short.pdf",plot_dir),width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

rfmix_pc_cor <- function(w=6,h=2){
	
	pcs <- fread(sprintf('%s/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt',meta_dir),sep='\t',header=T)
	global_ai <- fread(sprintf('%s/gtex_ai_global_all.txt',meta_dir),sep='\t',header=T)
	pcs[,SUBJID := sapply(FID, function(x) paste(unname(unlist(strsplit(x,'-')))[1:2],collapse='-'))]
	m <- merge(pcs, global_ai, by.x='SUBJID', by.y='ID')
	m[,UNK := NULL]
	m <- data.frame(m)
	rownames(m) <- m$SUBJID
	m$SUBJID <- NULL
	m$FID <- NULL
	m$IID <- NULL
	corr <- cor(m)
	g <- ggcorrplot(corr[1:5, 21:23],lab=T, outline.col='black') + 
		theme(panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.text.x = element_text(hjust=1,colour='black'),
			axis.text.y = element_text(colour='black')) +
		scale_x_discrete(labels=paste0('gPC',1:5))

	pdf(sprintf("%s/fig1b_rfmix_pc_cor.pdf",plot_dir),width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

local_block_trunc <- function(chr_start,chr_stop,w=7.5,h=2){

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
	dict <- data.frame(Chr=c(chr_start:chr_stop),Max=NA,Start=NA,Stop=NA)
	for (i in c(chr_start:chr_stop)){
		rows <- which(master$Chr == i)
		min <- min(master$Start[rows])
		master$Start[rows] <- master$Start[rows] - min # make each chromosome start at 0
		master$Stop[rows] <- master$Stop[rows] - min # shift all segments over by that amount
		max <- max(master$Stop[rows])
		dict$Max[dict$Chr==i] <- max
	}

	# figure out start and stop positions for each chromosome
	for (chr in c(chr_start:chr_stop)){
		if (chr == chr_start){
			dict$Start[dict$Chr==chr] <- 0
			dict$Stop[dict$Chr==chr] <- dict$Max[dict$Chr==chr]
		} else {
			dict$Start[dict$Chr==chr] <- dict$Stop[dict$Chr==chr-1] # make chromosomes start and stop on the same point
			dict$Stop[dict$Chr==chr] <- dict$Start[dict$Chr==chr] + dict$Max[dict$Chr==chr]
		}
	}

	# add start value to every "start" and "stop" in master
	for (chr in c(chr_start:chr_stop)){
		rows <- which(master$Chr == chr)
		master$Start[rows] <- master$Start[rows] + dict$Start[dict$Chr==chr]
		master$Stop[rows] <- master$Stop[rows] + dict$Start[dict$Chr==chr]
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

	master <- data.table(master)
	print(head(master))
	print(dict)

	print(chr_stop)
	if (chr_stop == 22){

		g <- ggplot(master[Chr %in% c(chr_start:chr_stop)]) + 
			geom_segment(aes(x=-Start, xend=-Stop, y=Y, yend=Y, colour=Pop), size=1, lineend="butt") +
			geom_hline(yintercept=yint, size=0.5) +
			geom_point(x=NA,y=NA,shape=22,aes(fill=Pop),colour='black',size=5) +
			theme_bw() +
			theme(panel.border = element_rect(size=1, colour='black'),
		          panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank(), 
		          axis.ticks.y = element_blank(),
		          axis.ticks.x = element_blank(),
		          axis.text.x = element_blank(),
		          #axis.text.y = element_text(size = 4, colour='black'),
		          axis.text.y = element_text(colour='black', hjust=1, vjust=0.5, size=10),
		          plot.title = element_blank(),
		          #plot.title = element_text(hjust = 0.5, size = 14, colour='black'),
		          legend.text = element_text(colour='black', size=10),
		          legend.title = element_blank(),
		          legend.position = 'bottom',
		          axis.title.y = element_blank(),
		          axis.title.x = element_blank(),
		          legend.margin=margin(t = -4.5, unit='mm')) +
			scale_colour_manual(name="Local Ancestry", 
				values = c(ASN="green3",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
				guide = 'none') +
			scale_fill_manual(name="Local Ancestry", 
				values = c(ASN="green3",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
				breaks=c('AFR','EUR','ASN',"UNK"),
				labels=c(ASN="East Asian (7.3%)",EUR="European (30.3%)",AFR="African (61.2%)",UNK="Unknown (<1%)")) +
			geom_vline(xintercept=-dict$Stop[1:nrow(dict)-1], size=0.5) + 
			scale_x_discrete(limits=-((dict$Start+dict$Stop)/2),
		        labels=c(chr_start:chr_stop),
		        breaks=-((dict$Start+dict$Stop)/2)) + 
			scale_y_continuous(breaks=breaks, labels=ids, expand = c(0,0), limits=c(0.5,234.5)) +
		  	guides(shape = guide_legend(override.aes = list(size = 1))) +
		  	coord_flip()

	} else {

		g <- ggplot(master[Chr %in% c(chr_start:chr_stop)]) + 
			geom_segment(aes(x=-Start, xend=-Stop, y=Y, yend=Y, colour=Pop), size=1, lineend="butt") +
			geom_hline(yintercept=yint, size=0.5) +
			geom_point(x=NA,y=NA,shape=22,aes(fill=Pop),colour='black',size=5) +
			theme_bw() +
			theme(panel.border = element_rect(size=1, colour='black'),
		          panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank(), 
		          axis.ticks.y = element_blank(),
		          axis.ticks.x = element_blank(),
		          axis.text.x = element_blank(),
		          #axis.text.y = element_text(size = 4, colour='black'),
		          axis.text.y = element_text(colour='black', hjust=1, vjust=0.5, size=10),
		          plot.title = element_blank(),
		          #plot.title = element_text(hjust = 0.5, size = 14, colour='black'),
		          legend.title = element_blank(),
		          legend.position = 'none',
		          axis.title.y = element_blank(),
		          axis.title.x = element_blank()) +
			scale_colour_manual(name="Local Ancestry", 
				values = c(ASN="green3",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
				guide = 'none') +
			scale_fill_manual(name="Local Ancestry", 
				values = c(ASN="green3",EUR="#0000FF",AFR="#FF9900",UNK="gray"), 
				guide = 'none') +
			geom_vline(xintercept=-dict$Stop[1:nrow(dict)-1], size=0.5) + 
			scale_x_discrete(limits=-((dict$Start+dict$Stop)/2),
		        labels=c(chr_start:chr_stop),
		        breaks=-((dict$Start+dict$Stop)/2)) + 
			scale_y_continuous(breaks=breaks, labels=ids, expand = c(0,0), limits=c(0.5,234.5)) +
		  	coord_flip()

	}

	pdf(sprintf("%s/block_haplo_trunc_%s_%s.pdf",plot_dir,chr_start,chr_stop),width=w,height=h)
	print(g)
	dev.off()

	return(g)

}
# g <- local_block_trunc(1,4)
# g <- local_block_trunc(19,22,h=1)

qq <- function(w=5,h=4){

	if (!file.exists(paste0(master_data_dir,'qq.RData'))){

		#load(paste0(master_data_dir,'egenes_master.RData'))
		load(paste0(master_data_dir,'allp_master.RData'))

		#egenes_master[,pval_nominal := as.numeric(pval_nominal)]

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
	print(head(qq))
	q <- ggplot(qq, aes(x=global, y=lava, colour=factor(tissue))) +
		geom_point(size=0.4,alpha=0.7) +
		scale_colour_manual(values=tissuecols, labels=tissuelabs, name='GTEx v8 tissue', guide='none') +
		labs(x=expression('GlobalAA'~italic(P)*'-value quantile'), y=expression('LocalAA'~italic(P)*'-value quantile')) +
		geom_abline(slope=1, intercept=0, linetype='dashed')  +
		theme_classic() +
		coord_cartesian(ylim=c(10,30),xlim=c(10,25))

	# too much info to save as PDF 
	jpeg(paste0(plot_dir,"/fig2a_qq.jpg"),res=300,quality=100,width=w,height=h,units='in')
	print(q)
	dev.off()

	return(q)

}

egene_discovery <- function(h=3,w=5,cutoff=1e-6){

	load(sprintf('%s/egenes_master-20200326.RData',master_data_dir))
	egenes_master <- egenes_master[pval_nominal_global < cutoff | pval_nominal_local < cutoff]

	counts <- data.table(TISSUE=tissues)
	counts[, total_global := 0]
	counts[, total_local := 0]

	for (t in tissues){
		
		sub <- egenes_master[tissue == t]

		t_local <- length(unique(sub[pval_nominal_local < cutoff, gene_id]))
		t_global <- length(unique(sub[pval_nominal_global < cutoff, gene_id]))

		counts[TISSUE == t, total_global := t_global]
		counts[TISSUE == t, total_local := t_local]

	}

	egene_discovery <- ggplot(counts, aes(x=total_global, y=total_local, fill=factor(TISSUE))) +
		geom_point(colour='black', shape=21, size=2) +
		geom_abline(slope=1, intercept=0, linetype='dashed') +
		scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue', limits=tissues) +
		labs(x='N GlobalAA eGenes', y='N LocalAA eGenes') +
		theme_classic() 

	pdf(paste0(plot_dir,"/fig2b_egene_discovery_",cutoff,".pdf"),width=w,height=h)
	print(egene_discovery)
	dev.off()

	return(egene_discovery)

}

venn <- function(h=3.5,w=5,cutoff=1e-6){

	load(sprintf('%s/egenes_master-20200326.RData',master_data_dir))
	egenes_master <- egenes_master[pval_nominal_local < cutoff | pval_nominal_global < cutoff]

	counts <- data.table(TISSUE=tissues)
	counts[, N_diff_lead_snp := 0]
	counts[, N_same_lead_snp := 0] 
	counts[, N_global_only := 0]
	counts[, N_local_only := 0]
	counts[, total_global := 0]
	counts[, total_local := 0]

	for (t in tissues){
		
		sub <- egenes_master[tissue == t]

		# unique to each method
		global_only <- length(unique(sub[pval_nominal_global < cutoff & pval_nominal_local > cutoff, gene_id]))
		local_only <- length(unique(sub[pval_nominal_local < cutoff & pval_nominal_global > cutoff, gene_id]))

		# total 
		t_local <- length(unique(sub[pval_nominal_local < cutoff, gene_id]))
		t_global <- length(unique(sub[pval_nominal_global < cutoff, gene_id]))

		# same eGene, diff lead variant 
		sub = sub[pval_nominal_global < cutoff & pval_nominal_local < cutoff]
		n_same_lead_snp = length(unique(sub[LD == 1 | overlapping_lead_variants == 1, gene_id]))
		n_diff_lead_snp = length(unique(sub[LD < 1 & overlapping_lead_variants == 0, gene_id]))

		counts[TISSUE == t, N_diff_lead_snp := n_diff_lead_snp]
		counts[TISSUE == t, N_same_lead_snp := n_same_lead_snp] 
		counts[TISSUE == t, N_global_only := global_only]
		counts[TISSUE == t, N_local_only := local_only]
		counts[TISSUE == t, total_global := t_global]
		counts[TISSUE == t, total_local := t_local]

	}

	dt <- counts[,.(TISSUE, N_diff_lead_snp, N_same_lead_snp, N_global_only, N_local_only)]

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
	                             labels=c(N_diff_lead_snp=' Same eGene,\n different lead eVariant',
	                                      N_same_lead_snp=' Same eGene,\n same lead eVariant',
	                                      N_global_only=' eGene in GlobalAA only ',
	                                      N_local_only=' eGene in LocalAA only'),
	                             breaks=c('N_local_only', 'N_global_only', 'N_same_lead_snp', 'N_diff_lead_snp')) +
	         theme_classic() +
	         scale_x_discrete(labels=shortlab,limits=tissues) +
	         theme(axis.title.x=element_blank(),
	               axis.text.x=element_text(colour='black'),
	               legend.title=element_blank(),
	               legend.position=c(0.5,0.87),
	               legend.key.height=unit(1.75, 'lines')) +
	         labs(y='N eGenes') +
	         ylim(c(0,max(dt[,local_y])+0.5*(max(dt[,local_y])))) +
	         guides(fill=guide_legend(ncol=2))

	pdf(paste0(plot_dir,"/fig2c_venn_",cutoff,".pdf"),width=w,height=h)
	print(g)
	dev.off()

	return(g)

}

pval_distn <- function(cutoff=1e-6,w=6,h=5){

	load(sprintf('%s/egenes_master-20200326.RData',master_data_dir))

	# ID genes unique to one method at given cutoff 
	la <- egenes_master[,.(gene_id, tissue, pval_nominal_global, pval_nominal_local, LD, overlapping_lead_variants)]	
	la <- la[ (pval_nominal_global < cutoff & pval_nominal_local > cutoff) | (pval_nominal_global > cutoff & pval_nominal_local < cutoff) ]
	la <- unique(la) 
	la <- la[order(overlapping_lead_variants, decreasing=T)]
	la2 <- la[abs(-log10(pval_nominal_global)-(-log10(pval_nominal_local)))>2]
	la2[,local_better := ifelse(pval_nominal_global > pval_nominal_local, 1, 0)]

	g <- ggplot(la, aes(x=-log10(pval_nominal_global), y=-log10(pval_nominal_local))) +
		geom_rect(ymin=-Inf, ymax=6, xmin=6, xmax=Inf, fill=methodcols[['global']]) +
		geom_rect(ymin=6, ymax=Inf, xmin=-Inf, xmax=6, fill=methodcols[['LAVA']]) +
		geom_vline(xintercept=6) +
		geom_hline(yintercept=6) +
		geom_abline(linetype='dashed') +
		geom_point(alpha=1,aes(shape=factor(overlapping_lead_variants), fill=tissue),colour='black', size=1) +
		theme_bw() +
		labs(x=expression('GlobalAA'~italic(P)*'-value (-log10)'), y=expression('LocalAA'~italic(P)*'-value (-log10)')) +
		scale_shape_manual(values=c('0'=21,'1'=24),name='Same lead\neVariant', labels=c('0'='No','1'='Yes')) +
		scale_fill_manual(values=tissuecols, guide='none') +
		theme(legend.position=c(0.8,0.81),
			panel.grid = element_blank()) +
		geom_abline(intercept=2, slope=1, linetype='dotted') +
		geom_abline(intercept=-2, slope=1, linetype='dotted') +
		guides(shape = guide_legend(override.aes = list(size=2))) +
		geom_rug(size=0.1,alpha=0.5) + 

	pdf(paste0(plot_dir,"/fig2d_uniq_egene_pval_distn_",cutoff,".pdf"),width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

gtex_r2 <- function(w=4,h=3){
	
	# plot_gtex_la_r2.R, r2_gtex_lead_snps.R
	#gtex_r2 <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/unique_lead_snp_r2.txt',sep='\t',header=FALSE)
	gtex_r2 = fread(sprintf('%s/unique_lead_snp_r2.txt',master_data_dir), sep='\t', header=T)

	g <- ggplot(gtex_r2, aes(V2)) +
		geom_histogram(bins=40,colour='black',fill='white') +
		theme_classic() + 
		labs(x='Variance in GTEx eVariant explained by LA',y='Frequency') +
		scale_x_continuous(limits=c(0,1), breaks=c(0.0,0.25,0.50,0.75,1.0), labels=c(0.0,0.25,0.50,0.75,1.0)) 

	pdf(sprintf('%s/fig4a_gtex-la-r2.pdf',plot_dir),width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

gtex_esnps_r2_high <- function(w=10,h=7){
	
	# d <- fread(cmd='zcat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/by-tissue/gtex_lead_snps_r2_la.merged.tsv.gz', sep='\t', header=T)
	# d[,tissue := sapply(tissue, function(x) gsub('_leadsnp.*','',x))]

	# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2-all_signif_pairs.RData')

	# m <- merge(d, gtex_esnp_r2[,.(tissue,variant_id,gene_id,pval_nominal)], by=c('tissue','variant_id'))
	# m <- unique(m)

	# msub <- m[local_rsq > 0.7]
	# m <- NULL
	# print(nrow(msub))
	# print(head(msub))

	# msub[,pval_nominal := -log10(pval_nominal)]

	# gtex_coloc <- fread(cmd='zcat /mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/gtex_coloc/results_p4_gt_0.5.txt.gz',sep='\t',header=F)
	# colnames(gtex_coloc) <- c('tissue','trait','gene_id','p4')
	# print(head(gtex_coloc))

	# m <- merge(msub, gtex_coloc, by=c('gene_id','tissue'), all.x=TRUE)

	# save(m, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/m_merged-gtex-eqtl-coloc.RData')
	load(sprintf('%s/m_merged-gtex-eqtl-coloc.RData',master_data_dir))
	
	m <- m[order(local_rsq, decreasing=T)]
	print(head(m))
	print(nrow(m))
	#write.table(m,'~/gtex-admix/metadata/gtex_esnps_high_r2_with_la.txt',sep='\t',col.names=T,row.names=F,quote=F)

	gene_names <- fread(sprintf('%s/gene_name_map.tsv',meta_dir), sep='\t', header=TRUE)
	m[,gene_stable_id := sapply(gene_id, function(x) gsub("\\..*","",x))]
	m <- merge(m, gene_names, by='gene_stable_id', all.m=T)

	labs <- c(Skin_Sun_Exposed_Lower_leg='SE skin',
		Pancreas='Pancreas',
		Adipose_Visceral_Omentum='Visc. adipose',
		Whole_Blood='Blood',
		Artery_Tibial='Tibial artery')

	# remove redundant labels
	m <- m[order(local_rsq, decreasing=T)]
	
	m[!is.na(trait),label := gsub('_',' ',paste0(gene_name,': ',trait))]
	m[duplicated(label),label := NA]
	print(m[!is.na(label)])

	print(all_labs[names(all_labs)%in%unique(m[!is.na(label),tissue])])
	print(unique(m[!is.na(label),tissue]))

	# shuffle
	m <- m[sample(1:nrow(m), nrow(m), replace=F)]

	g <- ggplot(data=m, aes(x=local_rsq, y=pval_nominal)) +
		geom_point(data=m,alpha=0.7,aes(colour=tissue)) +
		theme_classic() +
		labs(y=expression('GTEx eQTL'~italic(P)*'-value (-log10)'), x='Variance in GTEx eVariant explained by LA') +
		scale_colour_manual(values=tissuecols,breaks=unique(m[!is.na(label),tissue]),limits=unique(m[!is.na(label),tissue]),labels=all_labs[names(all_labs)%in%unique(m[!is.na(label),tissue])]) +
		
		geom_label_repel(data=m[!is.na(label)],
			aes(label=label,
				fill=tissue,
				colour=tissue),
			nudge_x=0.68-m[!is.na(label),local_rsq],
			#nudge_x=-250,
			direction='y',
			hjust=1,
			seed=1,
			alpha=0.4,
			size=3,
			show.legend=F) +
		
		geom_label_repel(data=m[!is.na(label)],
			aes(label=label,
				colour=tissue),
			nudge_x=0.68-m[!is.na(label),local_rsq],
			#nudge_y=-250,
			direction='y',
			hjust=1,
			seed=1,
			fill=NA,
			alpha=1,
			size=3,
			show.legend=F) +
		
		geom_label_repel(data=m[!is.na(label)],
			aes(label=label),
			nudge_x=0.68-m[!is.na(label),local_rsq],
			#nudge_y=-250,
			direction='y',
			hjust=1,
			seed=1,
			fill=NA,
			alpha=1,
			size=3,
			colour='black',
			label.size=NA,
			segment.size=NA,
			show.legend=F) + 

		theme(legend.position='bottom',
			legend.title=element_blank(),
			legend.text=element_text(size=11),
			axis.title=element_text(size=12),
			legend.margin=margin(t = -3, r=1, l=1, unit='mm'),
			axis.text.x = element_text(size=c(12,9,9,12), colour='black'),
			axis.text.y = element_text(colour='black')) +
		scale_fill_manual(values=tissuecols,guide='none') +
		guides(label = F, colour = guide_legend(override.aes=list(shape=19,size=3,alpha=0.7),ncol=3)) +
		scale_y_continuous(breaks=c(0,10,20,30,40,50)) +
		scale_x_continuous(limits=c(0.5,1), breaks=c(0.7,0.8,0.9,1.0)) +
		geom_vline(xintercept=0.7,linetype='dashed') 

	pdf(sprintf('%s/fig4b_gtex-high-la-r2-with-coloc-v2.pdf',plot_dir),width=w,height=h)
	print(g)
	dev.off()

	return(g)
}

fix_label = function(x){
	y = gsub('self_reported','SR',x)
	y = gsub('diagnosed_by_doctor','DBD',y)
	y = gsub('imputed_','',y)
	y = gsub('_et_al_2016','',y)
	y = gsub('.*_[0-9]','UKB',y)
	return(y)
}

coloc <- function(w=11,h=8){

	# output of compile_coloc.R --> now reformat_coloc.R
	# already excludes loci with lead variants in LD of 1
	load(sprintf('%s/master_coloc-1e-04-20200410.RData',master_data_dir))

	# master_coloc[,line_colour := ifelse(clpp_h4_global > clpp_h4_local, 'global', 'LAVA')]
	#master_coloc[, label:= ifelse((clpp_h4_global < 0.60 & clpp_h4_local > 0.75) | (clpp_h4_local < 0.60 & clpp_h4_global > 0.75), paste0(feature,':',gwas_trait), NA)]

	#master_coloc[,label := NULL]

	# # remove genes that don't have different lead SNPs
	# load("/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/merged/mx_merged-all-annotations-v2.RData")
	# la <- mx[,.(gene_id, tissue, lead_snp_D)]
	# master_coloc <- merge(master_coloc, la, by.x=c('feature','tissue'), by.y=c('gene_id','tissue'))
	# master_coloc <- master_coloc[lead_snp_D < 1]

	gene_names <- fread(sprintf('%s/gene_name_map.tsv',meta_dir), sep='\t', header=TRUE)

	master_coloc[, gene_stable_id := sapply(gene_id, function(x) unname(unlist(strsplit(x, '[.]')))[1])]
	master_coloc <- merge(master_coloc, gene_names, by='gene_stable_id', all.x=TRUE)

	master_coloc[,gwas_trait := sapply(gwas_trait, fix_label)]
	master_coloc[,label_global:= ifelse(clpp_h4_local < 0.50 & clpp_h4_global > 0.5 & (clpp_h4_global - clpp_h4_local) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]
	master_coloc[,label_local:= ifelse(clpp_h4_local > 0.5 & clpp_h4_global < 0.5 & (clpp_h4_local - clpp_h4_global) > 0.3, gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

	writeLines(unique(master_coloc[!is.na(label_local)|!is.na(label_global),gwas_trait]))
	#master_coloc[,label := ifelse(((clpp_h4_local < 0.60 & clpp_h4_global > 0.75) | (clpp_h4_local > 0.75 & clpp_h4_global < 0.60)), gsub('_',' ',paste0(gene_name,': ',gwas_trait)), NA)]

	g <- ggplot() +

		geom_rect(data=data.frame(xmin=0.5,xmax=1,ymin=0,ymax=0.5), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=0.3, fill=unname(methodcols['global'])) + # global is better
		geom_rect(data=data.frame(xmin=0,xmax=0.5,ymin=0.5,ymax=1), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=0.3, fill=unname(methodcols['LAVA'])) + # local is better

		geom_rect(data=data.frame(xmin=0.5,xmax=1,ymin=0,ymax=0.5), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), colour=unname(methodcols['global']),fill=NA) + # global is better
		geom_rect(data=data.frame(xmin=0,xmax=0.5,ymin=0.5,ymax=1), aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), colour=unname(methodcols['LAVA']),fill=NA) + # local is better
		
		geom_segment(data=data.frame(x=0,xend=1,y=0,yend=1),aes(x=x,y=y,xend=xend,yend=yend),colour='black',linetype='dashed') +
		geom_segment(data=data.frame(x=0,xend=1,y=0,yend=0),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
		geom_segment(data=data.frame(x=0,xend=1,y=1,yend=1),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
		geom_segment(data=data.frame(x=0,xend=0,y=1,yend=0),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
		geom_segment(data=data.frame(x=1,xend=1,y=0,yend=1),aes(x=x,y=y,xend=xend,yend=yend),size=0.5) +
	
		geom_point(data=master_coloc, aes(x=clpp_h4_global, y=clpp_h4_local, colour=tissue),alpha=0.7) +
		geom_point(data=master_coloc[!is.na(label_local) | !is.na(label_global)],alpha=1,aes(fill=tissue,clpp_h4_global, y=clpp_h4_local),colour='black',shape=21) +

		# local better, fill
		geom_label_repel(data=master_coloc[!is.na(label_local)],
			nudge_x=1.05 - master_coloc[!is.na(label_local),clpp_h4_global],
			direction="y",
			nudge_y=0.1,
			hjust=0,
			ylim=c(0.68,1),
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local,fill=tissue,colour=tissue),
			seed=1,
			alpha=0.4,
			size=3,
			show.legend=F) + 

		# local better, no fill
		geom_label_repel(data=master_coloc[!is.na(label_local)],
			nudge_x=1.05 - master_coloc[!is.na(label_local),clpp_h4_global],
			direction="y",
			nudge_y=0.1,
			hjust=0,
			ylim=c(0.68,1),
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local,colour=tissue),
			seed=1,
			fill=NA,
			size=3,
			show.legend=F) + 
		
		# local better, black text 
		geom_label_repel(data=master_coloc[!is.na(label_local)],
			nudge_x=1.05 - master_coloc[!is.na(label_local),clpp_h4_global],
			direction="y",
			nudge_y=0.1,
			hjust=0,
			ylim=c(0.68,1),
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_local),
			seed=1,
			fill=NA,
			colour='black',
			label.size=NA,
			segment.size=NA,
			size=3,
			show.legend=F) + 

		geom_label_repel(data=master_coloc[!is.na(label_global)],
			nudge_x=1.05 - master_coloc[!is.na(label_global),clpp_h4_global],
			nudge_y=-0.3,
			direction="y",
			hjust=0,
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global,fill=tissue,colour=tissue),
			seed=1,
			alpha=0.4,
			ylim=c(0,0.68),
			size=3,
			show.legend=F) + 

		geom_label_repel(data=master_coloc[!is.na(label_global)],
			nudge_x=1.05 - master_coloc[!is.na(label_global),clpp_h4_global],
			nudge_y=-0.3,
			direction="y",
			hjust=0,
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global,colour=tissue),
			seed=1,
			fill=NA,
			ylim=c(0,0.68),
			size=3,
			show.legend=F) + 

		geom_label_repel(data=master_coloc[!is.na(label_global)],
			nudge_x=1.05 - master_coloc[!is.na(label_global),clpp_h4_global],
			nudge_y=-0.3,
			direction="y",
			hjust=0,
			aes(x=clpp_h4_global,y=clpp_h4_local,label=label_global),
			seed=1,
			fill=NA,
			colour='black',
			label.size=NA,
			segment.size=NA,
			ylim=c(0,0.68),
			size=3,
			show.legend=F) + 

		scale_colour_manual(values=tissuecols, labels=tissuelabs, name='Tissue', limits=tissues) +
		scale_fill_manual(values=tissuecols, labels=tissuelabs, name='Tissue',guide='none') +

		labs(x='GlobalAA colocalization probability (PP4)', y='LocalAA colocalization probability (PP4)') +
		theme_classic() +

		scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1.85), expand = c(0,0)) +
		scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1), expand = c(0,0)) +
		guides(colour = guide_legend(override.aes=list(alpha=0.7,shape=19,size=3),ncol=6),fill=FALSE,label=F) + 
		theme(panel.grid=element_blank(),
			legend.position='top',
			legend.justification='left',
			legend.title=element_blank(),
			legend.margin=margin(b = -2, unit='mm'),
			legend.direction='horizontal',
			axis.title.x=element_text(hjust=0.15),
			legend.text=element_text(size=11)
			# axis.title=element_text(size=14)
			) +
		annotate("text",x=0.75,y=0.025,hjust=0.5,label='Stronger GlobalAA coloc.') +
		annotate("text",x=0.25,y=0.975,hjust=0.5,label='Stronger LocalAA coloc.') 
	
	png(paste0(plot_dir,'/fig4a_COLOC-local-global-p4-scatter-bytissue-labelled-v3.png'), width=w,height=h,unit='in',res=300)
	print(g)
	dev.off()

	return(g)
}

fst_distn <- function(w=3,h=5){
	# # compile_coloc.R
	# load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_all.RData')

	# # add eur-afr Fst
	# eur_afr <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/high_local_fst/EUR_AFR.weir.fst',sep='\t',header=T)
	# colnames(eur_afr) <- c('CHROM','POS','EUR_AFR_Fst')

	# eur_afr[EUR_AFR_Fst == '-nan', EUR_AFR_Fst := 0]
	# eur_afr = eur_afr[EUR_AFR_Fst != 'WEIR_AND_COCKERHAM_FST']

	# eur_afr[,EUR_AFR_Fst := as.numeric(EUR_AFR_Fst)]
	# eur_afr[,CHROM := as.character(CHROM)]
	# eur_afr[,POS := as.integer(POS)]
	# master[,CHROM := as.character(CHROM)]
	# master[,POS := as.integer(POS)]

	# master <- merge(master, eur_afr, by=c('CHROM','POS'))
	# master <- unique(master)

	# eur_max <- round(max(master[grepl('EUR',which_max_fst),max_fst]),digits=2)
	# afr_max <- round(max(master[grepl('AFR',which_max_fst),max_fst]),digits=2)
	# eur_afr_fst <- round(max(master[,EUR_AFR_Fst], na.omit=T),digits=2)

	# master <- melt(master, measure.vars=c('max_fst','EUR_AFR_Fst'))
	# print(head(master))

	# master[grepl('AFR',which_max_fst)&variable=='max_fst',group := 'Africa']
	# master[grepl('EUR',which_max_fst)&variable=='max_fst',group := 'Europe']
	# master[variable=='EUR_AFR_Fst',group := 'Between-continent Fst']

	# labels <- data.table(label=c(sprintf('Max within-EUR Fst: %s\nMax within-AFR Fst: %s',eur_max,afr_max),
	# 							sprintf('Max between-continent Fst\n(EUR vs AFR): %s',eur_afr_fst)),
	# 					variable=c('Max. within-continent Fst','Between-continent Fst'))
	# print(head(master))

	# master[variable == 'EUR_AFR_Fst', variable := 'Between-continent Fst']
	# master[variable == 'max_fst', variable := 'Max. within-continent Fst']

	# save(master, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_fst-distn-plot.RData')

	load(sprintf('%s/master_fst-distn-plot.RData',master_data_dir))
	print(max(master[group=='Europe'&variable=="Max. within-continent Fst",value], na.rm=T))
	print(max(master[group=='Africa'&variable=='Max. within-continent Fst',value], na.rm=T))

	master[variable == 'Max. within-continent Fst', name := 'F[ST*","*within]']
	master[variable == 'Between-continent Fst', name := 'F[ST*","*between]']
	# master[variable == "Max. within-continent Fst", title := 'Max. within-continent F[ST]']
	# master[variable == "Between-continent Fst", title := 'Between-continent F[ST]']

	print(master[value > 0.98])

	print(max(master[variable=='Between-continent Fst',value], na.rm=T))

	master <- master[complete.cases(master)]

	labels <- data.table(label='Max: 0.95',name='F[ST*","*between]')
	print(labels)

	print(head(master))

	g <- ggplot() +
		geom_histogram(data=master, aes(value, fill=group), bins=40, colour='black') +
		geom_label(data=labels, aes(label=label),x=Inf,y=Inf,hjust=1.2,vjust=1.4,label.size=NA) +
		theme_bw() + 
		labs(x=bquote('F'['ST']),y='Frequency') +
		scale_fill_manual(values=c(Europe="#0000FF",Africa="#FF9900",'Between-continent Fst'='white'),
			limits=c('Europe','Africa'),
			labels=c(Europe='Europe\n(max: 0.61)',Africa='Africa\n(max: 0.99)')) +
		facet_wrap(~name, ncol=1, labeller = label_parsed, scales='free_y') +
		scale_x_continuous(limits=c(0,0.5)) +
		theme(
			#strip.text = element_text(size=14),
			panel.border = element_rect(colour = "black", fill=NA, size=1),
			strip.background = element_rect(size=1, colour='black',fill=NA),
			legend.position=c(0.65,0.30),
			legend.title=element_blank(),
			#legend.spacing.y = unit(1.2, 'lines'),
			#axis.title=element_text(size=14),
			#legend.text=element_text(size=14),
			panel.grid = element_blank(),
			legend.key.width=unit(0.5,"line"),
			legend.key.heigh=unit(1.5,"line"))

	pdf(paste0(plot_dir,'/fst_distn-merged-by-superpop.pdf'), width=w, height=h)
	print(g)
	dev.off()

	return(g)
}

max_fst_pop <- function(w=5,h=4){
	# compile_coloc.R
	#load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_all.RData')
	load(sprintf('%s/master_all_fst.RData',master_data_dir))

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
		theme_classic() + 
		labs(x='Pairwise subpopulations',y=bquote('Frequency of F'['ST'*","*'within']~'(%)')) +
		scale_x_discrete(breaks=order_of_labels, limits=order_of_labels) +
		scale_fill_identity(guide='legend',name='Super-\npopulation',labels=c('#0000FF'='Europe','#FF9900'='Africa')) +
		theme(legend.position=c(0.65,0.2),
			axis.text.y=element_text(colour='black'),
			#axis.title=element_text(size=14),
			#legend.text=element_text(size=14),
			#legend.title=element_text(size=14)
			) +
		coord_flip() +
		scale_y_continuous(expand = c(0,0),limits=c(0,max(table[,N]+0.5)))

	pdf(paste0(plot_dir,'/fst_by_pair-all.pdf'), width=w, height=h)
	print(g)
	dev.off()

	return(g)
}

# need to redo this 

# regression_res <- function(w=3,h=4){

# 	dt = fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/merged-regression-results6.txt',sep='\t',header=T)
# 	print(dt)
# 	dt <- dt[variable!='(Intercept)']

# 	dt[,ymin := estimate - stderr]
# 	dt[,ymax := estimate + stderr]

# 	dt1 = dt

# 	# rename facets
# 	dt1[variable == 'max_within_continent_fst_local', name := 'F[ST*","*within]']
# 	dt1[variable == 'max_within_continent_fst_global', name := 'F[ST*","*within]']
# 	dt1[variable == 'WEIR_AND_COCKERHAM_FST_local', name := 'F[ST*","*between]']
# 	dt1[variable == 'WEIR_AND_COCKERHAM_FST_global', name := 'F[ST*","*between]']

# 	dt1[variable == 'max_within_continent_fst_local', colour := 'LAVA']
# 	dt1[variable == 'max_within_continent_fst_global', colour := 'global']
# 	dt1[variable == 'WEIR_AND_COCKERHAM_FST_local', colour := 'LAVA']
# 	dt1[variable == 'WEIR_AND_COCKERHAM_FST_global', colour := 'global']

# 	# # rename facets
# 	# dt1[variable == 'within', name := 'Delta*F[ST*","*within]']
# 	# dt1[variable == 'between', name := 'Delta*F[ST*","*between]']
# 	# # dt[variable == 'within_global', name := 'F[ST*","*within*","*GlobalAA]']
# 	# # dt[variable == 'between_local', name := 'F[ST*","*between*","*LocalAA]']
# 	# # dt[variable == 'within_local', name := 'F[ST*","*within*","*LocalAA]']

# 	dt1[pvalue < 0.01, label := '*']
# 	dt1[pvalue < 1e-3, label := '**']
# 	dt1[pvalue < 1e-5, label := '***']
# 	dt1[pvalue >= 0.05, label := '']

# 	# get_range <- function(var){
# 	#   sub = dt[variable == var]
# 	#   range = abs(max(sub[,ymax]) - min(sub[,ymin]))
# 	#   return(range)
# 	# }

# 	# for (var in unique(dt1[,variable])){
# 	#   print(var)
# 	#   print(get_range(var))
# 	#   dt1[variable==var, range := get_range(var)]
# 	# }

# 	dt1[,label_pos := ymax + 0.02]

# 	print(dt1)

# 	# change facet colors

# 	p <- ggplot(dt1, aes(x=colour)) +
# 	  geom_bar(aes(y=estimate, alpha=name, fill=colour),position='dodge',stat='identity',colour='black') +
# 	  theme_classic() + 
# 	  geom_errorbar(aes(ymin=ymin, ymax=ymax, group=name),width=0.1,position=position_dodge(width=0.9)) +
# 	  labs(x='Lead eVariant',y='Coefficient estimate') +
# 	  geom_text(aes(label=label,y=label_pos,group=name),show.legend=F,position=position_dodge(width=0.9),vjust=0.7) +
# 	  theme(axis.text=element_text(colour='black'),
# 			legend.title=element_blank(),
# 			legend.text=element_text(hjust=0),
# 			legend.position=c(0.75,0.15),
# 			axis.title.y=element_blank()) +
# 	  scale_fill_manual(values=methodcols, labels=methodlabs) +
# 	  scale_x_discrete(labels=methodlabs)+
# 	  guides(text=FALSE,fill=F) +
# 	  scale_alpha_manual(values=c('F[ST*","*within]'=0.6,'F[ST*","*between]'=1), 
# 	  	labels=c('F[ST*","*within]'=expression(F[ST*","*within]),'F[ST*","*between]'=expression(F[ST*","*between])),
# 	  	limits=c('F[ST*","*within]', 'F[ST*","*between]')) +
# 	  coord_flip()


# 	pdf(paste0(plot_dir,'/coloc-fst-regression-6.pdf'), width=w, height=h)
# 	print(p)
# 	dev.off()

# 	return(p)
# }

ve_expr_ancestry <- function(w=5,h=5){

	# variance_explained_ancestry_gene_expression.R; plot_var_expression_ancestry.R

	i <-1 
	dt_list <- list()
	for (t in tissues){
		#dt <- fread(sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/tss/%s_var_explained_by_ancestry.tsv',t),sep='\t',header=T)
		dt = fread(sprintf('%s/expr_ve/%s_var_explained_by_ancestry.tsv',master_data_dir, t),sep='\t',header=T)
		dt[,tissue := t]
		dt_list[[i]] <- dt
		i <- i + 1
	}

	merged <- rbindlist(dt_list)

	print(head(merged))
	print(t.test(merged[,variance_explained_ga], merged[,variance_explained_la]))

	# shuffle rows
	merged <- merged[sample(1:nrow(merged),nrow(merged),replace=F)]

	# add gene names
	gene_names <- fread(sprintf('%s/gene_name_map.tsv',meta_dir), sep='\t', header=TRUE)
	merged <- merge(merged, gene_names, by.x='gene', by.y='gene_stable_id', all.x=TRUE)

	print(mean(merged[,variance_explained_la]))
	print(sd(merged[,variance_explained_la]))
	print(merged[order(variance_explained_la,decreasing=T)][1:20])
	print(merged[gene_name=='TBC1D3L'])
	print(t.test(merged[gene_name=='TBC1D3L',variance_explained_la],merged[gene_name=='TBC1D3L',variance_explained_ga]))

	collapsed <- merged[,list(mean_la=mean(variance_explained_la),mean_ga=mean(variance_explained_ga)),by=c('gene','gene_name')]
	collapsed[,diff := mean_la - mean_ga]
	collapsed <- collapsed[order(diff,decreasing=T)]
	print(collapsed[1:20])

	merged[,colour := tissue]
	merged[gene_name=='TBC1D3L', colour := 'black']
	merged[,size := ifelse(gene_name=='TBC1D3L', 2, 1)]
	merged[,alpha := ifelse(gene_name=='TBC1D3L', 2, 1)]

	m1 = merged[gene_name=='TBC1D3L']
	m2 = merged[gene_name != 'TBC1D3L']
	merged = rbind(m2, m1)

	g <- ggplot(merged, aes(x=variance_explained_ga*100, y=variance_explained_la*100)) +
		geom_point(aes(colour=factor(colour), 
			size=factor(size), 
			fill=factor(tissue), 
			alpha=factor(alpha)), 
		shape=21) +
		theme_classic() + 
		labs(x='% of residual expr. VE by GA',y='% of residual expr. VE by LA') +
		scale_fill_manual(values=tissuecols, guide='none') +
		scale_colour_manual(values=c(tissuecols, black='black'),
			breaks='black',
			labels='TBC1D3L') +
		theme(legend.title = element_blank(), 
			legend.position=c(0.83,0.67),
			legend.background=element_blank()) +
		geom_abline(linetype='dashed') +
		scale_size_manual(values=c('1'=1,'2'=2), guide='none') +
		scale_alpha_manual(values=c('1'=0.5,'2'=1), guide='none') +
		guides(colour = guide_legend(override.aes=list(shape=21,size=2)))

	png(paste0(plot_dir,'/var-explained.png'),width=w,height=h,unit='in',res=300)
	print(g)
	dev.off()

	return(g)
}

########################################################################################################################################
# supplementary figures

admix_per_tissue <- function(w=4,h=6){
	
	aa <- sprintf('%s/aa-per-tissue.tsv',meta_dir)
	aa <- fread(aa, sep='\t', header=TRUE)

	aa[,colour := ifelse(tissue %in% tissues, 'black', NA)]

	print(head(aa))

	aa <- aa[order(N_AA,decreasing=T)]
	limits <- rev(aa[,tissue])

	# change tissue colors
	tcols <- tissuecols
	tcols[!names(tcols) %in% tissues] <- 'gray'

	g <- ggplot(aa,aes(x=tissue,y=N_AA,fill=tissue)) +
	geom_bar(stat="identity",aes(colour=colour)) +
	geom_text(
		aes(label = frac, y = N_AA + 2),
		#position = position_dodge(0.9),
		hjust = 0,
		vjust = 0.5,
		size = 2.5) +
	scale_fill_manual(values=tcols,guide='none') +
	scale_colour_manual(values=c(black='black','NA'=NA),guide='none') +
	scale_x_discrete(labels=labs,limits=limits) +
	labs(y="N admixed individuals (117AX)") +
	theme_classic() +
	theme(axis.text.y = element_text(hjust=1,vjust=0.5,
		# size=tissize[match(limits, names(tissize))],
		# face=face[match(limits, names(face))],
		colour='black'),
		axis.title.y=element_blank(),
		legend.margin=margin(t = 0, unit='mm')) +
    scale_y_continuous(expand = c(0,0), limits=c(0,115)) +
    coord_flip()

	pdf(sprintf("%s/admix_per_tissue.pdf",supp_dir),width=w,height=h)
	print(g)
	dev.off()
}

tg_gtex_pca <- function(){

	#infile <- '/mnt/lab_data/montgomery/nicolerg/pca-subset/gtexV8-PC12-df.tsv'
	infile = sprintf('%s/gtexV8-PC12-df.tsv',master_data_dir)
	dat <- fread(infile, sep='\t', header=TRUE)

	group.colors <- c(ESN="#FFCD00",
		GWD="#FFB900",
		LWK="#CC9933",
		MSL="#E1B919",
		YRI="#FFB933",
		ACB="#FF9900",
		ASW="#FF6600",
		CLM="#CC3333",
		MXL="#E10033",
		PEL="#FF0000",
		PUR="#CC3300",
		CDX="#339900",
		CHB="#ADCD00",
		CHS="#00FF00",
		JPT="#008B00",
		KHV="#00CC33",
		CEU="#0000FF",
		FIN="#00C5CD",
		GBR="#00EBFF",
		IBS="#6495ED",
		TSI="#00008B",
		BEB="#8B008B",
		GIH="#9400D3",
		ITU="#B03060",
		PJL="#E11289",
		STU="#FF00FF")

	breaks <- c("ESN",
		"GWD",
		"LWK",
		"MSL",
		"YRI",
		"ACB",
		"ASW",
		"CLM",
		"MXL",
		"PEL",
		"PUR",
		"CDX",
		"CHB",
		"CHS",
		"JPT",
		"KHV",
		"CEU",
		"FIN",
		"GBR",
		"IBS",
		"TSI",
		"BEB",
		"GIH",
		"ITU",
		"PJL",
		"STU")

	g <- ggplot(data=dat) +
		geom_point(data=subset(dat, !pop=="GTEx"), aes(x=PC1, y=PC2, colour=pop, shape=dataset), alpha=1, size=1) + 
		theme_classic() + 
		theme(legend.position="right",
			plot.title=element_text(hjust = 0.5),
			legend.text=element_text(size=12),
			legend.title=element_text(size=15),
			axis.text=element_text(size=10),
			legend.key.size=unit(18,"point"),
			axis.title=element_text(size=15)) +
		scale_colour_manual(values=group.colors,name="1000 Genomes",breaks=breaks) +
		geom_point(data=subset(dat, pop=="GTEx"), aes(x=PC1, y=PC2, shape=dataset), colour="black", size=1.5) +
		scale_shape_manual(values=c(20,2),name="Consortium") +
		xlab("gPC1") + ylab("gPC2") 

	pdf(sprintf('%s/gtex-tg-geno-PC12.pdf',supp_dir), width=8, height=6)
	print(g)
	dev.off()

}

# from locus_plots.R
# this needs to be rerun first 

# locus_plots <- function(){
	
# 	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master-locus_plots_v3.RData')
# 	#master[,p4_lab := sapply(p4_lab, function(x) gsub('p4','PP4',p4_lab))]
# 	print(head(master))

# 	# first just where Local is better (7 out of 23)
# 	local <- master[global_p4 < local_p4]
# 	g <- ggplot(local) +
# 		geom_point(data=local[!is.na(variant_id)&global_p4 < local_p4],colour='black',size=2,aes(x=pvalue_gwas,y=pvalue_eqtl)) +
# 		geom_point(aes(colour=method,x=pvalue_gwas,y=pvalue_eqtl),alpha=0.7,size=0.5) +
# 		scale_colour_manual(name='eQTL ancestry adjustment',values=c(Global=methodcols[['global']],Local=methodcols[['LAVA']]), labels=c(Global='GlobalAA',Local='LocalAA')) +
# 		theme_bw() +
# 		guides(color = guide_legend(override.aes = (list(size=2)), title.position = "top", 
# 									title.hjust = 0.5)) +
# 		theme(legend.position='top',
# 			legend.margin=margin(b = -2, unit='mm'),
# 			strip.text = element_text(size=8),
# 			panel.border = element_rect(colour = "black", fill=NA, size=1),
# 			strip.background = element_rect(size=1, colour='black',fill=NA),
# 			panel.grid = element_blank()) +
# 		labs(x=expression('GWAS'~italic(P)*'-value (-log10)'),y=expression('eQTL'~italic(P)*'-value (-log10)')) +
# 		facet_wrap(~description,scales='free',ncol=4,labeller = label_wrap_gen()) 
# 	g1 <- g + geom_text(data = unique(local,by=c('description','p4_lab')), mapping = aes(x = -Inf, y = Inf, label = p4_lab),hjust = -0.1,vjust = 1.2,size=2.5)

# 	png(sprintf('%s/coloc-locus-compare-wrapped-better-local.png',supp_dir),width=8,height=5,units='in',res=300)
# 	print(g1)
# 	dev.off()

# 	# now just where Global is better (16 out of 23)
# 	global <- master[global_p4 > local_p4]
# 	g <- ggplot(global) +
# 		geom_point(data=global[!is.na(variant_id)&global_p4 > local_p4],colour='black',size=2,aes(x=pvalue_gwas,y=pvalue_eqtl)) +
# 		geom_point(aes(colour=method,x=pvalue_gwas,y=pvalue_eqtl),alpha=0.7,size=0.5) +
# 		scale_colour_manual(name='eQTL ancestry adjustment',values=c(Global=methodcols[['global']],Local=methodcols[['LAVA']]), labels=c(Global='GlobalAA',Local='LocalAA')) +
# 		theme_bw() +
# 		guides(color = guide_legend(override.aes = (list(size=2)), title.position = "top", 
# 									title.hjust = 0.5)) +
# 		theme(legend.position='top',
# 			legend.margin=margin(b = -2, unit='mm'),
# 			strip.text = element_text(size=8),
# 			panel.border = element_rect(colour = "black", fill=NA, size=1),
# 			strip.background = element_rect(size=1, colour='black',fill=NA),
# 			panel.grid = element_blank()) +
# 		labs(x=expression('GWAS'~italic(P)*'-value (-log10)'),y=expression('eQTL'~italic(P)*'-value (-log10)')) +
# 		facet_wrap(~description,scales='free',ncol=4,labeller = label_wrap_gen()) 
# 	g2 <- g + geom_text(data = unique(global,by=c('description','p4_lab')), mapping = aes(x = -Inf, y = Inf, label = p4_lab),hjust = -0.1,vjust = 1.2, size=2.5)

# 	png(sprintf('%s/coloc-locus-compare-wrapped-better-global.png',supp_dir),width=8,height=9,units='in',res=300)
# 	print(g2)
# 	dev.off()

# 	# split into two pages (one where local better; one where global better)

# 	# pdf(sprintf("%s/coloc-locus-compare.pdf",supp_dir),width=15, height=15)
# 	# grid.arrange(g1, g2,
# 	# 	layout_matrix = rbind(c(1,1),
# 	# 						c(2,2),
# 	# 						c(2,2),
# 	# 						c(2,2))
# 	# )
# 	# dev.off()

# }

########################################################################################################################################
# put them together 

g_legend <- function(a.gplot){
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}

figure1 <- function(){

	g <- local_block_trunc(1,2)

	pdf(sprintf("%s/figure1-tmp.pdf",plot_dir),width=8, height=8)
	grid.arrange(genotype_pc(), admix_per_tissue_short(), rfmix_pc_cor(), g, ve_expr_ancestry(),
		layout_matrix = rbind(c(1,2),
								c(1,2),
								c(1,2),
								c(4,4),
								c(3,5),
								c(3,5))
	)
	dev.off()

	
	g <- local_block_trunc(21,22,h=1)
	# add 19-22 manually 
}

figure4 <- function(){

	gtex <- gtex_r2()
	high <- gtex_esnps_r2_high()
	#manually put them together 

	# pdf(sprintf("%s/figure4.pdf",plot_dir),width=8, height=4)
	# grid.arrange(gtex, high,
	# 	layout_matrix = rbind(c(1,1,2,2,2),
	# 							c(1,1,2,2,2))
	# )
	# dev.off()

}

figure2 <- function(){
	
	r <- rectGrob(gp=gpar(fill="white",colour='white'))

	# pdf(sprintf("%s/figure2-tmp.pdf",plot_dir),width=7.5,height=6)
	# grid.arrange(egene_discovery(cutoff=1e-06), r, venn(cutoff=1e-06), pval_distn(cutoff=1e-06),
	# 	layout_matrix = rbind(c(2,2,2,1,1,1,1),
	# 						c(2,2,2,1,1,1,1),
	# 						c(3,3,3,3,4,4,4),
	# 						c(3,3,3,3,4,4,4),
	# 						c(3,3,3,3,4,4,4))
	# )
	# dev.off()

	png(sprintf("%s/figure2.png",plot_dir),width=7.5,height=6,units='in',res=300)
	grid.arrange(egene_discovery(cutoff=1e-06), qq(), venn(cutoff=1e-06), pval_distn(cutoff=1e-06),
		layout_matrix = rbind(c(2,2,2,1,1,1,1),
							c(2,2,2,1,1,1,1),
							c(3,3,3,3,4,4,4),
							c(3,3,3,3,4,4,4),
							c(3,3,3,3,4,4,4))
	)
	dev.off()

}

figure3 <- function(cutoff=1e-6){

	coloc <- coloc()
	pop_fst <- fst_distn()
	max_pop <- max_fst_pop()
	reg <- regression_res()

	png(sprintf("%s/figure3.png",plot_dir),res=300, units='in', width=8.5, height=10)
	grid.arrange(coloc, pop_fst, max_pop, reg,
		layout_matrix=rbind(c(1,1,1,1,1,1,1,1,1),
							c(1,1,1,1,1,1,1,1,1),
							c(1,1,1,1,1,1,1,1,1),
							c(2,2,2,3,3,3,4,4,4),
							c(2,2,2,3,3,3,4,4,4))
	)
	dev.off()

}

# supplement
# tg_gtex_pca()
# admix_per_tissue()
# locus_plots()

args <- parser$parse_args()

c <- as.numeric(args$cutoff)

if(args$genotype_pc){
	g <- genotype_pc()
}
if(args$admix_per_tissue_short){
	g <- admix_per_tissue_short()
}
if(args$rfmix_pc_cor){
	g <- rfmix_pc_cor()
}
if(args$local_block_trunc){
	g <- local_block_trunc(1,4)
	g <- local_block_trunc(19,22,h=1)
}
if(args$qq){
	g <- qq()
}
if(args$egene_discovery){
	g <- egene_discovery(cutoff=c)
}
if(args$venn){
	g <- venn(cutoff=c)
}
if(args$pval_distn){
	g <- pval_distn(cutoff=c)
}
if(args$gtex_r2){
	g <- gtex_r2()
}
if(args$gtex_esnps_r2_high){
	g <- gtex_esnps_r2_high()
}
if(args$coloc){
	g <- coloc()
}
if(args$fst_distn){
	g <- fst_distn()
}
if(args$max_fst_pop){
	g <- max_fst_pop()
}
if(args$regression_res){
	g <- regression_res()
}
if(args$ve_expr_ancestry){
	g <- ve_expr_ancestry()
}
if(args$tg_gtex_pca){
	g <- tg_gtex_pca()
}
if(args$admix_per_tissue){
	g <- admix_per_tissue()
}
if(args$locus_plots){
	g <- locus_plots()
}
if(args$figure1){
	g <- figure1()
}
if(args$figure2){
	g <- figure2()
}
if(args$figure3){
	g <- figure3(cutoff=c)
}
if(args$figure4){
	g <- figure4()
}


