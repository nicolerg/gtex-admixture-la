# regulatory enrichment:
# 1. Take shared eGenes between LAVA and GLOBAL.
#	Look at subset where eSNP differs. 
#	Test enrichment of annotations.
# 2. Take top SNP for GEUVADIS 5% eGenes (EUR-> AFR and AFR->EUR) for global and LAVA.
#	Doesn’t matter if gene is significant or not. 
#	Look at where eSNP differs. 
#	Test enrichment of annotations.

indir="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/maf-filt/reg-enrich"

# eqtldir="/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/AFR-into-EUR"
# Rscript ~/gtex-admix/scripts/LAVA/1kg-eqtl/pval-correct-fastqtl.R ${eqtldir} "AFR.geu.permuted.txt.gz"
# Rscript ~/gtex-admix/scripts/LAVA/1kg-eqtl/pval-correct-fastqtl.R ${eqtldir} "EUR.geu.permuted.txt.gz"
# eqtldir="/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/EUR-into-AFR"
# Rscript ~/gtex-admix/scripts/LAVA/1kg-eqtl/pval-correct-fastqtl.R ${eqtldir} "AFR.geu.permuted.txt.gz"
# Rscript ~/gtex-admix/scripts/LAVA/1kg-eqtl/pval-correct-fastqtl.R ${eqtldir} "EUR.geu.permuted.txt.gz"

# Rscript ~/gtex-admix/scripts/geuvadis/reg-gene-sets.R

base="/mnt/lab_data/montgomery/nicolerg/local-eqtl/lung/aa-only/maf-filt/"

# # geneset = same genes, different lead SNP
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py ${base}/global.gtex.cutoff.eqtl.genes.tsv ${base}/lava.gtex.cutoff.eqtl.genes.tsv ${base}/gtex.diff.lead.esnp.txt "gtex.diff.lead.snp" &
# # geneset = top 100 genes from GEUVADIS 
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py ${base}/gtex.admix.global.egenes.txt.gz ${base}/gtex.admix.lava.egenes.txt.gz ${base}/top100.afr2eur.txt "top100.afr2eur" &
# # geneset = top 100 genes from GEUVADIS 
# python ~/gtex-admix/scripts/geuvadis/parse-regulatory-enrichment.py ${base}/gtex.admix.global.egenes.txt.gz ${base}/gtex.admix.lava.egenes.txt.gz ${base}/top100.eur2afr.txt "top100.eur2afr" &

# wait

wait
Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${base}/reg-enrich/global.regenrich.gtex.diff.lead.snp.tsv.gz ${base}/reg-enrich/lava.regenrich.gtex.diff.lead.snp.tsv.gz "gtex.diff.leadsnp"
Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${base}/reg-enrich/global.regenrich.top100.afr2eur.tsv.gz ${base}/reg-enrich/lava.regenrich.top100.afr2eur.tsv.gz "top100.afr2eur"
Rscript ~/gtex-admix/scripts/geuvadis/plot-reg-enrichment.R ${base}/reg-enrich/global.regenrich.top100.eur2afr.tsv.gz ${base}/reg-enrich/lava.regenrich.top100.eur2afr.tsv.gz "top100.eur2afr"
