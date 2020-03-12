library(ggplot2)
library(data.table)

# system("echo 'Reading files'")
# localfile <- "/mnt/lab_data/montgomery/nicolerg/admix/wholeblood-chr21subset-allpairs-localcov-hg19-500-noPCs.tsv"
# localcov <- "/mnt/lab_data/montgomery/nicolerg/admix/chr21_subsets/chr21.hg19.localcov.tsv"

# output <- fread(localfile, sep='\t', header=TRUE)
# cov <- read.delim(localcov, sep='\t', header=TRUE, stringsAsFactors=FALSE)
# glob <- read.delim('/users/nicolerg/gtex-admix/scripts/global-ids.txt',sep='\t',header=FALSE)

# colnames(glob) <- c("SUBJID","RACE")
# white <- nrow(glob[glob$RACE==3,])

# local_anc_frac <- data.frame(matrix(nrow=length(unique(cov$SNP_ID)),ncol=3,NA))
# colnames(local_anc_frac) <- c("SNP_ID","AFR_FRAC","ASN_FRAC")

# head(local_anc_frac)

# output$AFR_FRAC <- NA
# output$ASN_FRAC <- NA

# system("echo 'Calculating local ancestry proportions'")
# # find fraction of local ancestry for each snp
# for(i in 1:length(unique(cov$SNP_ID))){
# 	snp <- unique(cov$SNP_ID)[i]
# 	print(paste("Working on ",snp,sep=''))
# 	sub <- cov[grep(snp, cov$SNP_ID),4:7]
# 	afrsum <- sum(sub$AFR)
# 	asnsum <- sum(sub$ASN)
# 	eursum <- sum(sub$EUR) + 2*white
# 	unksum <- sum(sub$UNK)
# 	total <- sum(afrsum, asnsum, eursum, unksum)
# 	afr <- afrsum/total
# 	asn <- asnsum/total
# 	local_anc_frac$SNP_ID[i] <- snp
# 	local_anc_frac$AFR_FRAC[i] <- afr
# 	local_anc_frac$ASN_FRAC[i] <- asn
# 	output$ASN_FRAC[grep(snp,output$variant_id)] <- asn
# 	output$AFR_FRAC[grep(snp,output$variant_id)] <- afr
# 	print(head(local_anc_frac))
# }

# write.table(output, "/users/nicolerg/gtex-admix/output_updated.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

# output$afr_pval_nominal <- -log(output$afr_pval_nominal)
# output$asn_pval_nominal <- -log(output$asn_pval_nominal)

# write.table(local_anc_frac, "/users/nicolerg/gtex-admix/local_ancestry_fractions.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
# write.table(output, "/users/nicolerg/gtex-admix/output_updated.tsv", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

output <- read.delim("/users/nicolerg/gtex-admix/output_updated.tsv", sep='\t', header=TRUE)

system("echo 'Plotting'")

jpeg("/users/nicolerg/gtex-admix/local_frac_vs_pval_african.jpg", width=1600, height=1400, res=200)
ggplot(output, aes(x=afr_pval_nominal,y=AFR_FRAC)) + 
	geom_point() + 
	labs(x="Nominal P-Value of AFR Local Covariate (-log)",y="Fraction of Subjects Assigned AFR (per SNP)") +
	ylim(0.1,0.15)
dev.off()

jpeg("/users/nicolerg/gtex-admix/local_frac_vs_pval_asian.jpg", width=1600, height=1400, res=200)
ggplot(output, aes(x=asn_pval_nominal,y=ASN_FRAC)) + 
	geom_point() + 
	labs(x="Nominal P-Value of ASN Local Covariate (-log)",y="Fraction of Subjects Assigned ASN (per SNP)") +
	ylim(0.0114,0.0153)
dev.off()