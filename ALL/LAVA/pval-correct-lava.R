# DEPRECATED - use extract-egene-lava.py instead 


library(data.table)
library(qvalue)
library(ggplot2)

# args <- commandArgs(trailingOnly = TRUE)
# indir <- args[1]

#indir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/'
indir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/by-chr/maf-filt/'
global <- paste0(indir,"gtex.admix.global.egenes.txt.gz")
lava <- paste0(indir,"gtex.admix.lava.egenes.txt.gz")

lava <- fread(input=sprintf("zcat %s", lava), sep='\t', header=TRUE)
global <- fread(input=sprintf("zcat %s", global), sep='\t', header=TRUE)

correct <- function(df){
	df[,pval_bf := pval_nominal * count]
	df[pval_bf > 1, pval_bf := 1]
	df[,qval := qvalue(pval_bf)$qvalues]
	df <- df[qval < 0.05]
	return(df)
}

print(paste('N LAVA genes:',nrow(lava)))
print(paste('N global genes:',nrow(global)))

lava1 <- correct(lava)
print(paste('LAVA BF+BH:', nrow(lava1)))
global1 <- correct(global)
print(paste('Global BF+BH:', nrow(global1)))

# compare to GTEx v8 cutoffs

egenes <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/Lung.v8.egenes.txt.gz"
lung <- fread(input=sprintf("zcat %s",egenes), sep='\t', header=TRUE)

lungmap <- lung[,list(gene_id, pval_nominal_threshold)]

n <- rep(0,4)

lavam <- merge(lava,lungmap,by="gene_id")
lavam <- lavam[pval_nominal < pval_nominal_threshold]
n[1] <- nrow(lavam)
print(paste('LAVA GTEx cutoff:', nrow(lavam)))
write.table(lavam,paste0(indir,'/lava.gtex.cutoff.eqtl.genes.tsv'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)

globalm <- merge(global,lungmap,by="gene_id")
globalm <- globalm[pval_nominal < pval_nominal_threshold]
n[2] <- nrow(globalm)
print(paste('Global GTEx cutoff:', nrow(globalm)))
write.table(lavam,paste0(indir,'/global.gtex.cutoff.eqtl.genes.tsv'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)

sharedg <- merge(globalm,lavam,by="gene_id")
print(paste('Shared egenes:',nrow(sharedg)))
n[3] <- nrow(sharedg)
print(paste('Shared esnp:',nrow(sharedg[variant_id.x==variant_id.y])))
n[4] <- nrow(sharedg[variant_id.x==variant_id.y])

# df <- data.frame(group=c("LAVA","global"),
# 	all=c(6864,6701),
# 	egene=c(3940,3940),
# 	esnp=c(2340))

# pdf("~/gtex-admix/plots/figures/eqtl-discoveries.pdf",height=4,width=5)
# ggplot(df, aes(x=group)) +
# 	geom_bar(aes(y=all,fill="red"),stat="identity") +
# 	geom_bar(aes(y=egene,fill="blue"),stat="identity") +
# 	geom_bar(aes(y=esnp,fill="yellow"),stat="identity") +
# 	scale_fill_manual(name='',values=c("red"="lightblue","blue"="plum","yellow"="palegreen"),
# 		labels=c("red"="N total eGenes","blue"="N shared eGenes","yellow"="N shared pairs"),
# 		breaks=c("red","blue","yellow"),
# 		guide='legend') +
# 	labs(y="N eQTL genes") +
# 	theme_bw() +
# 	theme(axis.text.x=element_text(colour="black",size=12),
# 		axis.title.x=element_blank(),
# 		axis.title.y=element_text(size=12),
# 		legend.text=element_text(size=12)) +
# 	scale_x_discrete(limits=c("LAVA","global"))
# dev.off()

# df <- data.frame(x=c("LAVA","global","eGene","eSNP"),
# 	y=c(6864,6701,3940,2340))

df <- data.frame(x=c("LAVA","global","eGene","eSNP"),
	y=n)

df

pdf("~/gtex-admix/plots/figures/eqtl-discoveries-split.pdf",height=4,width=5)
ggplot(df, aes(x=x,y=y), colour="#99FF00") +
	geom_bar(stat="identity") +
	labs(y="N eQTL genes") +
	theme_bw() +
	theme(axis.text.x=element_text(colour="black",size=12),
		axis.title.x=element_blank(),
		axis.title.y=element_text(size=12)) +
	scale_x_discrete(limits=c("LAVA","global","eGene","eSNP"),
		labels=c(LAVA="LAVA",global="global",eGene="shared\neGenes",eSNP="shared\neGene-eSNP\npairs")) 
dev.off()

