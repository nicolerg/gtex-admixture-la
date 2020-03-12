library(ggcorrplot)
library(data.table)
library(plyr)

tissues <- c('Adipose_Subcutaneous',
	'Nerve_Tibial',
	'Artery_Tibial',
	'Muscle_Skeletal',
	'Lung',
	'Skin_Not_Sun_Exposed_Suprapubic')

for (t in tissues){
	
	gtex_covariates <- fread(sprintf('/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/%s.v8.covariates.txt',t),header=TRUE, sep='\t')
	gtex_peers <- gtex_covariates[grepl("InferredCov", `ID`)]

	new_peers <- fread(sprintf('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/%s/peer/%s.PEER_covariates.txt',t,t), header=TRUE, sep='\t')

	gtex_peers[,`ID` := paste0('gtex_v8_',`ID`)]
	new_peers[,`ID` := paste0('admixed_',`ID`)]

	gtex_peers <- gtex_peers[,colnames(gtex_peers)[colnames(gtex_peers) %in% colnames(new_peers)], with=FALSE]

	gtex_peers <- as.data.frame(gtex_peers)
	new_peers <- as.data.frame(new_peers)

	rownames(gtex_peers) <- gtex_peers$`ID`
	gtex_peers$`ID` <- NULL
	rownames(new_peers) <- new_peers$`ID`
	new_peers$`ID` <- NULL

	gtex_peers <- gtex_peers[,match(colnames(gtex_peers), colnames(new_peers))]

	peers <- rbind(gtex_peers, new_peers)
	print(head(peers))

	subj <- colnames(peers)
	peer <- rownames(peers)

	peers <- data.frame(t(peers))
	colnames(peers) <- peer
	rownames(peers) <- subj

	corr <- cor(peers)

	corr <- corr[1:20, 61:75]
	print(head(corr))

	g <- ggcorrplot(corr)
	g <- g + labs(title=t)

	pdf(sprintf('~/gtex-admix/plots/new_peers/%s_peer_corr.pdf',t),height=6, width=9)
	print(g)
	dev.off()

}