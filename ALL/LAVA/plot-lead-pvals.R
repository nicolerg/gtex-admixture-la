library(data.table)
library(ggplot2)

both <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Muscle_Skeletal/merged.eQTL.diff.leadSNP.pvalues.tsv"
lead <- "/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/Muscle_Skeletal/merged.eQTL.diff.leadSNP.tsv"

both <- fread(both, sep='\t', header=TRUE)
lead <- fread(lead, sep='\t', header=TRUE)

lead[,diff := abs(pval_nominal_global - pval_nominal_lava)]
both[,diff := ifelse(method_where_lead == "global", global_pval_nominal - lava_pval_nominal, lava_pval_nominal - global_pval_nominal)]

# order by diff
lead <- lead[order(diff)]
both <- both[order(diff)]
lead[,index:=1:nrow(lead)]
both[,index:=1:nrow(both)]

jpeg("~/gtex-admix/plots/muscle-leadsnp-pval-diff.jpg",width=1200,height=600,res=200)
ggplot(lead, aes(x=index, y=diff)) +
	geom_bar(stat="identity") +
	theme(axis.text.x=element_blank(),
		axis.title.x=element_blank()) +
	labs(y="abs(p-value difference between lead SNPs)", title="eGenes with different lead SNPs")
dev.off()

jpeg("~/gtex-admix/plots/muscle-leadsnp-match-pval-diff.jpg",width=1400,height=600,res=200)
ggplot(both, aes(x=index, y=diff, colour=factor(method_where_lead))) +
	geom_bar(stat="identity") +
	theme(axis.text.x=element_blank(),
		axis.title.x=element_blank()) +
	labs(y="p-value diff (where lead - where not lead)", title="Difference in p-values for same SNP between methods") +
	scale_colour_manual(values=c(lava="red",global="blue"),name="Method where lead")
dev.off()

jpeg("~/gtex-admix/plots/muscle-leadsnp-pval-diff-hist.jpg",width=1200,height=600,res=200)
ggplot(lead, aes(x=diff)) +
	geom_histogram(bins=100) +
	labs(x="abs(p-value difference between lead SNPs)", y="frequency", title="eGenes with different lead SNPs")
dev.off()

jpeg("~/gtex-admix/plots/muscle-leadsnp-match-pval-diff-hist.jpg",width=1200,height=600,res=200)
ggplot(both, aes(x=diff,fill=method_where_lead)) +
	geom_histogram(bins=100,position="identity",alpha=0.5) +
	labs(x="abs(p-value difference between lead SNPs)", y="frequency", title="eGenes with different lead SNPs") +
	scale_fill_manual(values=c(lava="red",global="blue"),name="Method where lead") 
dev.off()
