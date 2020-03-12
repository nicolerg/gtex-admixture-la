# compare effect sizes from permuted FastQTL in different populations

library(data.table)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]

indir <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/"
outdir <- indir

afr <- paste0(indir,"AFR.chr",chr,".merged.txt.gz")
eur <- paste0(indir,"EUR.chr",chr,".merged.txt.gz")

afr <- fread(input=sprintf("zcat %s", afr), sep=' ', header=FALSE)
eur <- fread(input=sprintf("zcat %s", eur), sep=' ', header=FALSE)

cols <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
colnames(afr) <- paste("AFR",cols,sep='.')
colnames(eur) <- paste("EUR",cols,sep='.')
afr$AFR.bonferroni = p.adjust(afr$AFR.bpval, method="bonferroni")
eur$EUR.bonferroni = p.adjust(eur$EUR.bpval, method="bonferroni")

# 1. ID of the tested molecular phenotype (in this particular case, the gene ID)
# 2. Number of variants tested in cis for this phenotype
# 3. MLE of the shape1 parameter of the Beta distribution
# 4. MLE of the shape2 parameter of the Beta distribution
# 5. Dummy [To be described later]
# 6. ID of the best variant found for this molecular phenotypes (i.e. with the smallest p-value)
# 7. Distance between the molecular phenotype - variant pair
# 8. The nominal p-value of association that quantifies how significant from 0, the regression coefficient is
# 9. The slope associated with the nominal p-value of association [only in version > v2-184]
# 10. A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype.
# 11. A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis.

merged <- merge(afr, eur, by.x="AFR.pid", by.y="EUR.pid")
colnames(merged)[colnames(merged)=="AFR.pid"] <- "pid"

paste0("N genes: ",(nrow(merged)))
paste0("N matching var: ",length(grep("TRUE",merged$EUR.sid == merged$AFR.sid)))
paste0("N eQTL genes in both: ",nrow(merged[merged$EUR.bpval < 0.05 & merged$AFR.bpval < 0.05]))
paste0("N eQTL genes in AFR: ",nrow(merged[merged$AFR.bpval < 0.05]))
paste0("N eQTL genes in EUR: ",nrow(merged[merged$EUR.bpval < 0.05]))

# filter out eGenes with bonferroni pval == 1.0
merged <- merged[merged$EUR.bpval < 0.05 & merged$AFR.bpval < 0.05]
# merged[,c("EUR.slope","AFR.slope","EUR.bonferroni","AFR.bonferroni","EUR.bpval","AFR.bpval")]

# filter to include only eGenes with significantly different effect sizes 
mean <- mean(abs(c(merged$EUR.slope, merged$AFR.slope)))
var <- var(c(merged$EUR.slope, merged$AFR.slope))

merged <- merged[abs(merged$EUR.slope - merged$AFR.slope) >= mean]
merged[,c("EUR.slope","AFR.slope","EUR.bpval","AFR.bpval")]

write.table(merged, paste0(outdir,"eqtl.curated.geneset.chr",chr,".tsv"), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

# compare the curated gene set to GTEx eQTLs in each tissue
f <- "/mnt/lab_data/montgomery/nicolerg/geuvadis/eqtl/gtex-eqtl-genes-bytissue.tsv.gz"
mat <- fread(input=sprintf("zcat %s",f), sep='\t', header=TRUE)

merged[,pid:=gsub("[.].*","",pid)]
mat[,gene_id:=gsub("[.].*","",gene_id)]

m <- na.omit(match(merged[,pid], mat[,gene_id]))
matfilt <- mat[m]
matfilt[,gene_id:=NULL]
df <- data.frame(sort(colSums(matfilt)/nrow(merged), decreasing=TRUE))
colnames(df) <- "frac_eqtl"
df$rank <- 1:nrow(df)
df
