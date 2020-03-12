# plot correlation between nominal p-values with global covariates and nominal p-values of halplo, afr, and asn covariates from the regression with local covariates

library(ggplot2)

egenes_global <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/chr21_subsets/wholeblood.egenes.chr21.tsv"
egenes_local <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/wholeblood-local-chr21-egenes.tsv"

locale <- read.delim(egenes_local, sep='\t', header=TRUE, stringsAsFactors=FALSE)
globale <- read.delim(egenes_global, sep='\t', header=TRUE, stringsAsFactors=FALSE)

locale$global_pval <- -log10(globale$pval_nominal)
locale$asn_pval_nominal <- -log10(locale$asn_pval_nominal)
locale$afr_pval_nominal <- -log10(locale$afr_pval_nominal)
locale$pval_nominal <- -log10(locale$pval_nominal)

locale <- locale[locale$global_pval > 5,]

jpeg("/users/nicolerg/gtex-admix/egene-pval-correlation.jpg", width=1000, height=600, res=150)
ggplot(locale, aes(x=global_pval)) + theme_bw() + 
	geom_point(aes(y=pval_nominal, color="geno"),size=1.2,alpha=0.7) +
	geom_point(aes(y=asn_pval_nominal, colour="asn"),size=1.2,alpha=0.7) +
	geom_point(aes(y=afr_pval_nominal, colour="afr"),size=1.2,alpha=0.7) +
	scale_color_manual(breaks=c("geno", "asn", "afr"), labels = c("Genotype Dosage", "Asian", "African"), values=c("#FF9900", "#ADFF33", "gray"), name="Covariate") +
	labs(x="Nominal P-Value with Global Covariates (-log10)", y="Nominal P-Values with Local Ancestry (-log10)") 
dev.off()