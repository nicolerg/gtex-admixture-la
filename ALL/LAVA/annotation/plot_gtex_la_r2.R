library(data.table)
library(ggplot2)
library(ggrepel)

gtex_r2 <- fread('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/unique_lead_snp_r2.txt',sep='\t',header=FALSE)

# pdf('~/gtex-admix/plots/new_peers/gtex-la-r2.pdf',width=6,height=3)
# ggplot(gtex_r2, aes(V2)) +
# 	geom_histogram(bins=500) +
# 	theme_bw() + 
# 	labs(x='Variance in GTEx v8 eSNP genotypes explained by local ancestry',y='Frequency')
# dev.off()

colnames(gtex_r2) <- c('variant_id','r2')

if (!file.exists('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2-all_signif_pairs.RData')){

	egene_df <- list()
	i <- 1
	indir='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL'
	for (file in list.files(path=indir, pattern="signif_variant_gene_pairs")){
		egenes <- fread(cmd=sprintf('zcat %s/%s',indir,file))
		egenes <- egenes[,.(gene_id, variant_id, pval_nominal)]
		merged <- merge(egenes, gtex_r2, by='variant_id')
		t <- gsub('\\.v8.*','',file)
		print(t)
		merged[,tissue := t]

		egene_df[[i]] <- merged
		i <- i + 1
	}

	gtex_esnp_r2 <- data.table(rbindlist(egene_df))
	save(gtex_esnp_r2,file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2-all_signif_pairs.RData')

} else {
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2-all_signif_pairs.RData')
}

# if (!file.exists('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2.RData')){

# 	colnames(gtex_r2) <- c('variant_id','r2')

# 	egene_df <- list()
# 	i <- 1
# 	indir='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL'
# 	for (file in list.files(path=indir, pattern="v8\\.egenes\\.txt\\.gz")){
# 		egenes <- fread(cmd=sprintf('zcat %s/%s',indir,file))
# 		egenes <- egenes[pval_nominal < pval_nominal_threshold]
# 		egenes <- egenes[,.(gene_id, variant_id, pval_nominal)]
# 		merged <- merge(egenes, gtex_r2, by='variant_id')
# 		tissue <- gsub('\\.v8.*','',file)
# 		print(tissue)
# 		merged[,tissue := tissue]

# 		egene_df[[i]] <- merged
# 		i <- i + 1
# 	}

# 	gtex_esnp_r2 <- data.table(rbindlist(egene_df))
# 	save(gtex_esnp_r2,file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2.RData')


# } else {
# 	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_esnp_r2.RData')
# }

gtex_esnp_r2[,pval_nominal := -log10(pval_nominal)]

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

gtex_esnp_r2 <- gtex_esnp_r2[sample(1:nrow(gtex_esnp_r2))]

# jpeg('~/gtex-admix/plots/new_peers/gtex-la-r2-by-pval-2.jpg',width=1000,height=1000,res=200)
# ggplot(gtex_esnp_r2, aes(colour=tissue, x=pval_nominal, y=r2)) +
# 	geom_point(alpha=0.5) +
# 	theme_bw() +
# 	labs(x='GTEx eQTL p-value (-log10)', y='Variance in lead SNP genotype explained by LA') +
# 	theme(legend.position='none') +
# 	scale_colour_manual(values=tissuecols)
# dev.off()

# jpeg('~/gtex-admix/plots/new_peers/gtex-la-r2-by-pval-zoom-2.jpg',width=800,height=800,res=200)
# ggplot(gtex_esnp_r2[r2>0.70], aes(colour=tissue, x=pval_nominal, y=r2)) +
# 	geom_point(alpha=0.7) +
# 	theme_bw() +
# 	labs(x='GTEx eQTL p-value (-log10)', y='Variance in lead SNP genotype explained by LA') +
# 	scale_colour_manual(values=tissuecols) +
# 	theme(legend.position='none') 
# dev.off()

x <- gtex_esnp_r2[r2>0.7]
head(x)

# see if any of these SNPs have high colocalization (from GTEx GWAS subgroup results)

if(!file.exists('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_coloc-r2_all_esnps.RData')){
	coloc <- list()
	i <- 1
	for (file in list.files(path='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/_old_peer/annotation/gtex_gwas_coloc/results_enloc_priors')){
		c <- fread(cmd=sprintf('zcat %s/%s','/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/_old_peer/annotation/gtex_gwas_coloc/results_enloc_priors',file),sep='\t',header=TRUE)
		tis <- gsub(".*PM__","",file)
		t <- gsub("\\.txt\\.gz","",tis)
		gwas_trait <- gsub("__PM__.*","",file)
		print(gwas_trait)
		print(t)
		sub <- x[tissue==t]
		if(nrow(sub) > 0){
			merged <- merge(sub,c,by='gene_id')
			merged <- merged[p4>0.5]
			if(nrow(merged)>0){
				print(merged)
				merged[,trait:=gwas_trait]
				coloc[[i]] <- merged
				i <- i+1
			}
		}
		trait <- NULL
		tissue <- NULL
		sub <- NULL
		merged <- NULL 
	}

	gtex_coloc <- data.table(rbindlist(coloc))
	save(gtex_coloc, file='/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_coloc-r2_all_esnps.RData')
} else {
	load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/r2/gtex/results/gtex_coloc-r2_all_esnps.RData')
}

gtex_coloc <- gtex_coloc[p4>0.7]

m <- merge(x, gtex_coloc, by=c('variant_id','gene_id','pval_nominal','r2','tissue'), all.x=TRUE)
head(m)
nrow(m)
m[!is.na(p4)]

gene_names <- fread('~/gtex-admix/metadata/gene_name_map.tsv', sep='\t', header=TRUE)

m[,gene_stable_id := sapply(gene_id, function(x) gsub("\\..*","",x))]
m <- merge(m, gene_names, by='gene_stable_id')

labs <- c(Skin_Sun_Exposed_Lower_leg='Skin (sun-exposed)',
	Pancreas='Pancreas',
	Adipose_Visceral_Omentum='Visc. adipose',
	Whole_Blood='Whole blood',
	Artery_Tibial='Tibial artery')

m[order(r2,decreasing=TRUE)][1:20]

jpeg('~/gtex-admix/plots/new_peers/gtex-la-r2-by-pval-zoom.jpg',width=1200,height=800,res=200)
ggplot(data=m, aes(x=pval_nominal, y=r2)) +
	geom_point(data=m,alpha=0.7,aes(colour=tissue)) +
	theme_bw() +
	labs(x='GTEx eQTL p-value (-log10)', y='Variance in lead SNP genotype explained by LA') +
	scale_colour_manual(values=tissuecols,breaks=unique(m[!is.na(trait),tissue]),limits=unique(m[!is.na(trait),tissue]),labels=labs) +
	geom_label_repel(data=m[!is.na(trait)],
		aes(label=gsub('_',' ',paste0(gene_name,': ',trait)),
			fill=tissue),
		nudge_x=40-m[!is.na(trait),pval_nominal],
		nudge_y=0.15,
		direction='y',
		hjust=0,
		seed=1,
		alpha=0.4,
		size=3,
		ylim=c(0.76,NA)) +
	geom_label_repel(data=m[!is.na(trait)],
		aes(label=gsub('_',' ',paste0(gene_name,': ',trait))),
		nudge_x=40-m[!is.na(trait),pval_nominal],
		nudge_y=0.15,
		direction='y',
		hjust=0,
		seed=1,
		fill=NA,
		alpha=1,
		size=3,
		ylim=c(0.76,NA)) +
	xlim(c(0,120)) +
	theme(legend.position='top',
		legend.title=element_blank(),
		legend.text=element_text(size=10),
		legend.margin=margin(b = -3, unit='mm'),) +
	guides(colour = guide_legend(override.aes=list(shape=15,size=4),ncol=5),fill=FALSE) +
	scale_fill_manual(values=tissuecols)
dev.off()


