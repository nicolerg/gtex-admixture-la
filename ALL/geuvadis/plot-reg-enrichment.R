library(ggplot2)
library(data.table)
library(reshape)

args <- commandArgs(trailingOnly=TRUE)
global <- args[1]
local <- args[2]
suffix <- args[3]

globalreg <- fread(input=sprintf("zcat %s", global),sep='\t',header=TRUE)
lavareg <- fread(input=sprintf("zcat %s", local),sep='\t',header=TRUE)

print(suffix)
paste0('N global eQTLs: ',nrow(globalreg))
paste0('N LAVA eQTLs: ',nrow(lavareg))

elements <- c(
	"enhancer",
	"enhancer_activity",
	"promoter",
	"promoter_activity",
	"promoter_flanking_region",
	"promoter_flanking_region_activity",
	"open_chromatin_region",
	"CTCF_binding_site",
	"TF_binding_site")

comp <- data.table(element=elements,LAVA=0,global=0)

for(el in elements){
	if(grepl("activity",el)){
		n_lava = nrow(lavareg[get(el) == "ACTIVE"])
		n_global = nrow(globalreg[get(el) == "ACTIVE"])
	} else {
		n_lava = nrow(lavareg[get(el) != "NA"])
		n_global = nrow(globalreg[get(el) != "NA"])
	}
	comp[ element == el, LAVA := n_lava ]
	comp[ element == el, global := n_global ]
}

n_lava <- sum(comp[,LAVA])
n_global <- sum(comp[,global])

paste('LAVA el:',sum(comp[,LAVA]))
paste('global el:',sum(comp[,global]))

# comp[,LAVA := LAVA/nrow(lavareg)]
# comp[,global := global/nrow(globalreg)]

# comp[,LAVA := LAVA/nrow(lavareg)]
# comp[,global := global/nrow(globalreg)]
comp <- melt(comp, id="element")

labs <- c(
	"enhancer" = "enhancer",
	"enhancer_activity" = "enhancer\nactivity",
	"promoter" = "promoter",
	"promoter_activity" = "promoter\nactivity",
	"promoter_flanking_region" = "promoter\nflanking\nregion",
	"promoter_flanking_region_activity" = "promoter\nflanking\nregion\nactivity",
	"open_chromatin_region" = "open\nchromatin\nregion",
	"CTCF_binding_site" = "CTCF\nbinding\nsite",
	"TF_binding_site" = "TF\nbinding\nsite")

pdf(paste0("~/gtex-admix/plots/lava.global.regenrich.",suffix,".pdf"),width=8,height=4)
ggplot(comp, aes(x=factor(element),y=value,fill=factor(variable))) +
	geom_bar(stat="identity",position="dodge") +
	theme_bw() +
	labs(y="N annotations") +
	scale_fill_manual(name="Method",
		values=c(LAVA="firebrick",global="steelblue"),
		labels=c(LAVA=paste0('LAVA (N = ',n_lava,')'),
			global=paste0('Global (N = ',n_global,')'))) +
	scale_x_discrete(labels = labs) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x=element_blank())
dev.off()

globalreg[,colour:="global"]
lavareg[,colour:="lava"]
master <- rbind(globalreg,lavareg)

# # tss_distance
# jpeg(paste0("~/gtex-admix/plots/lava.global.tssdist.",suffix,".jpg"),res=200,width=800,height=600)
# ggplot(data=master) +
# 	geom_density(aes(x=tssdist,colour=factor(colour),fill=factor(colour)),alpha=0.5) +
# 	theme_bw() +
# 	scale_colour_manual(name="Method",values=c(global="steelblue",lava="firebrick")) +
# 	scale_fill_manual(name="Method",values=c(global="steelblue",lava="firebrick")) +
# 	labs(x="TSS distance", y="Density")
# dev.off()

# allp <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/chr1/lung-LAVA-global-merged-allpairs-chr1.tsv.gz"
# globalannot <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/reg-enrich/global-egene-regulatory-enrichment.tsv.gz"
# lavalannot <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/reg-enrich/LAVA-egene-regulatory-enrichment.tsv.gz"

# globalreg <- fread(input=sprintf("zcat %s", globalannot),sep='\t',header=TRUE)
# lavareg <- fread(input=sprintf("zcat %s", lavalannot),sep='\t',header=TRUE)
# allp <- fread(input=sprintf("zcat %s", allp),sep='\t',header=TRUE)

# cutoff <- 0.05/(nrow(allp))

# mergedglobal <- merge(globalreg, allp, by.x=c("gene","variant"), by.y=c("gene_id","variant_id"))
# mergedlocal <- merge(lavareg, allp, by.x=c("gene","variant"), by.y=c("gene_id","variant_id"))
# mergedglobal <- mergedglobal[global_pval_nominal < cutoff]
# mergedlocal <- mergedlocal[lava_pval_nominal < cutoff]

# elements <- c("open_chromatin_region",
# 	"enhancer_activity",
# 	"promoter_activity",
# 	"enhancer",
# 	"promoter",
# 	"promoter_flanking_region",
# 	"CTCF_binding_site",
# 	"promoter_flanking_region_activity",
# 	"TF_binding_site")

# comp <- data.table(element=elements,LAVA=0,global=0)

# for(el in elements){
# 	if(grepl("activity",el)){
# 		n_lava = nrow(mergedlocal[get(el) == "ACTIVE"])
# 		n_global = nrow(mergedglobal[get(el) == "ACTIVE"])
# 	} else {
# 		n_lava = nrow(mergedlocal[get(el) != "NA"])
# 		n_global = nrow(mergedglobal[get(el) != "NA"])
# 	}
# 	comp[ element == el, LAVA := n_lava ]
# 	comp[ element == el, global := n_global ]
# }

# comp
# comp

# comp[, LAVA := LAVA / nrow(mergedlocal) ]
# comp[, global := global / nrow(mergedglobal) ]

# # plot fraction in global vs LAVA for each element

# comp <- melt(comp, id="element")

# labs <- c(
# 	"enhancer" = "enhancer",
# 	"enhancer_activity" = "enhancer\nactivity",
# 	"promoter" = "promoter",
# 	"promoter_activity" = "promoter\nactivity",
# 	"promoter_flanking_region" = "promoter\nflanking\nregion",
# 	"promoter_flanking_region_activity" = "promoter\nflanking\nregion\nactivity",
# 	"open_chromatin_region" = "open\nchromatin\nregion",
# 	"CTCF_binding_site" = "CTCF\nbinding\nsite",
# 	"TF_binding_site" = "TF\nbinding\nsite")

# pdf("~/gtex-admix/plots/lava-global-reg-enrichment-allgenes.pdf",width=8,height=4)
# ggplot(comp, aes(x=factor(element),y=value,fill=factor(variable))) +
# 	geom_bar(stat="identity",position="dodge") +
# 	theme_bw() +
# 	labs(y="Fraction of shared eQTL genes (N = 72)") +
# 	scale_fill_manual(name="Method",values=c(LAVA="firebrick",global="steelblue")) +
# 	scale_x_discrete(labels = labs) +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x=element_blank())
# dev.off()
# # # tss_distance
# # jpeg("~/gtex-admix/plots/lava-global-tss-dist.jpg",res=200,width=800,height=600)
# # ggplot(data=NULL,aes(x=tss_distance)) +
# # 	geom_density(data=mergedglobal,aes(colour="Global",fill="Global"),alpha=0.5) +
# # 	geom_density(data=mergedlocal,aes(colour="LAVA",fill="LAVA"),alpha=0.5) +
# # 	theme_bw() +
# # 	scale_colour_manual(name="Method",values=c(Global="steelblue",LAVA="firebrick")) +
# # 	scale_fill_manual(name="Method",values=c(Global="steelblue",LAVA="firebrick")) +
# # 	labs(x="TSS distance", y="Density")
# # dev.off()

# # now do it for overlapping genes only ----------------------------------------------------------------------
# globalannot <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/reg-enrich/global-egene-regulatory-enrichment.tsv.gz"
# lavalannot <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/reg-enrich/LAVA-egene-regulatory-enrichment.tsv.gz"
# lavae <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/chr1/gtex.admix.lava.egenes.chr1.txt.gz"
# globale <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/chr1/gtex.admix.global.egenes.chr1.txt.gz"

# # gene_id count   pval_nominal    tss_distance    slope   variant_id

# globalreg <- fread(input=sprintf("zcat %s", globalannot),sep='\t',header=TRUE)
# lavareg <- fread(input=sprintf("zcat %s", lavalannot),sep='\t',header=TRUE)
# ge <- fread(input=sprintf("zcat %s", globale),sep='\t',header=TRUE)
# le <- fread(input=sprintf("zcat %s", lavae),sep='\t',header=TRUE)

# gmerged <- merge(globalreg,ge,by.x=c("gene","variant"),by.y=c("gene_id","variant_id"))
# print(paste0("N global egenes: ",nrow(gmerged)))
# lmerged <- merge(lavareg,le,by.x=c("gene","variant"),by.y=c("gene_id","variant_id"))
# print(paste0("N LAVA egenes: ",nrow(lmerged)))

# master <- merge(gmerged, lmerged, by="gene", suffixes=c(".global",".lava"))
# print(paste0("N overlapping egenes: ",nrow(master)))

# elements <- c(
# 	"enhancer",
# 	"enhancer_activity",
# 	"promoter",
# 	"promoter_activity",
# 	"promoter_flanking_region",
# 	"promoter_flanking_region_activity",
# 	"open_chromatin_region",
# 	"CTCF_binding_site",
# 	"TF_binding_site")

# comp <- data.table(element=elements,LAVA=0,global=0)

# for(el in elements){
# 	if(grepl("activity",el)){
# 		n_lava = nrow(master[get(paste0(el,".lava")) == "ACTIVE"])
# 		n_global = nrow(master[get(paste0(el,".global")) == "ACTIVE"])
# 	} else {
# 		n_lava = nrow(master[get(paste0(el,".lava")) != "NA"])
# 		n_global = nrow(master[get(paste0(el,".global")) != "NA"])
# 	}
# 	comp[ element == el, LAVA := n_lava ]
# 	comp[ element == el, global := n_global ]
# }

# comp
# comp

# comp[, LAVA := LAVA / nrow(master) ]
# comp[, global := global / nrow(master) ]

# # plot fraction in global vs LAVA for each element

# comp <- melt(comp, id="element")

# comp

# labs <- c(
# 	"enhancer" = "enhancer",
# 	"enhancer_activity" = "enhancer\nactivity",
# 	"promoter" = "promoter",
# 	"promoter_activity" = "promoter\nactivity",
# 	"promoter_flanking_region" = "promoter\nflanking\nregion",
# 	"promoter_flanking_region_activity" = "promoter\nflanking\nregion\nactivity",
# 	"open_chromatin_region" = "open\nchromatin\nregion",
# 	"CTCF_binding_site" = "CTCF\nbinding\nsite",
# 	"TF_binding_site" = "TF\nbinding\nsite")

# pdf("~/gtex-admix/plots/lava-global-reg-enrichment-samegenes.pdf",width=8,height=4)
# ggplot(comp, aes(x=factor(element),y=value,fill=factor(variable))) +
# 	geom_bar(stat="identity",position="dodge") +
# 	theme_bw() +
# 	labs(y="Fraction of shared eQTL genes (N = 72)") +
# 	scale_fill_manual(name="Method",values=c(LAVA="firebrick",global="steelblue")) +
# 	scale_x_discrete(labels = labs) +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x=element_blank())
# dev.off()
# # tss_distance
# jpeg("~/gtex-admix/plots/lava-global-tss-dist-samegenes.jpg",res=200,width=800,height=600)
# ggplot(data=master) +
# 	geom_density(aes(x=tss_distance.global,colour="Global",fill="Global"),alpha=0.5) +
# 	geom_density(aes(x=tss_distance.lava,colour="LAVA",fill="LAVA"),alpha=0.5) +
# 	theme_bw() +
# 	scale_colour_manual(name="Method",values=c(Global="steelblue",LAVA="firebrick")) +
# 	scale_fill_manual(name="Method",values=c(Global="steelblue",LAVA="firebrick")) +
# 	labs(x="TSS distance", y="Density")
# dev.off()
