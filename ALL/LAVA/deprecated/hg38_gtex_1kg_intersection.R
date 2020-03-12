library(data.table)

# DEPRECATED. use format_run_admixture for vcftools/bcftools implementation

# Get intersection of variants (not LD-pruned or subset in any way, except for chr1 and chr3)

args <- commandArgs(trailingOnly=TRUE)
chrom <- args[1]

outdir <- '/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixture/merged_vcf_hg38'

kg <- sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixture/rand_1_3/tg_hg38_chr%s_rand.vcf.gz',chrom)
gtex <- sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixture/gtex_hg38/rand/gtex_hg38_phased_maf01_rand_chr%s.vcf',chrom)
# kg <- sprintf('/mnt/lab_data/montgomery/shared/1KG/hg38/ALL.chr%s_GRCh38.genotypes.20170504.vcf.gz',chrom)
# gtex <- sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/geno/hg38/gtex.admixed.MAC10.phased.chr%s.hg38.vcf.gz',chrom)

# read in files 
kg <- fread(cmd=sprintf('zcat %s',kg), skip='#CHROM', sep='\t', header=T)
gtex <- fread(gtex, sep='\t', header=T)

# retain only admixed individuals in GTEx
gtex_admixed <- fread('/mnt/lab_data/montgomery/nicolerg/admixed-ancestry/gtex-admixed0.9.txt',sep='\t',header=F)
keep <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', gtex_admixed[,V1])
print(keep)
gtex <- gtex[,keep,with=F]

# prepend 'chr' to 1kg VCF
kg[,`#CHROM` := paste0('chr',`#CHROM`)]

# remove columns we don't want
kg[ ,c("QUAL","FILTER","INFO","FORMAT") := NULL]
gtex[ ,c("QUAL","FILTER","INFO","FORMAT") := NULL]

# merge
merged <- merge(kg, gtex, by=c('#CHROM','POS','REF','ALT'), suffixes=c('_kg','_gtex'))
# make a custom ID (allow mapping between GTEx ID and 1kg ID)
merged[, "ID" := paste0(ID_gtex, '-', ID_kg)]
merged[,c('ID_kg','ID_gtex') := NULL]

# add back in some columns 
filler <- c('QUAL','FILTER','INFO','FORMAT')
merged[, (filler) := NA]
ordered <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
setcolorder(merged, c(ordered, colnames(merged)[!colnames(merged) %in% ordered]))

print(sprintf('%s shared variants on chromosome %s', nrow(merged), chrom))

# write it out 
write.table(merged, sprintf('%s/gtex-1kg-hg38-merged-chr%s.vcf',outdir, chrom), sep='\t', col.names=T, row.names=F, quote=F)

