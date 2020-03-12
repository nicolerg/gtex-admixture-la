library(data.table)

ids <- '/mnt/lab_data/montgomery/nicolerg/rfmix-inputs/gtex-ids'
iddf <- read.table(ids, sep='\t', header=FALSE)
ids <- unname(unlist(iddf[,1]))
idlist <- rep(NA, 2*length(ids))
i <- 1
for(id in ids){
	idlist[i] <- paste0(id,'_A')
	idlist[i+1] <- paste0(id,'_B')
	i <- i + 2
}

dtlist <- vector("list", 22)
j <- 1
for(chr in 1:22){
	vit <- paste0('/mnt/lab_data/montgomery/nicolerg/rfmix-outputs/whole-output/GTEX-all-chr',chr,'.0.Viterbi.txt.gz')
	snp <- paste0('/mnt/lab_data/montgomery/nicolerg/rfmix-inputs/snp-bed-chr',chr,'_cM.bed')
	viterbi <- fread(input = sprintf("zcat %s",vit), sep=' ', header=FALSE)
	colnames(viterbi) <- idlist
	snpmap <- fread(snp, sep='\t', header=FALSE)
	snpmap <- snpmap[,1:4]
	colnames(snpmap) <- c("chr","start","end","variant_id")
	
	dt <- data.table(cbind(snpmap, viterbi))
	obj <- paste0('dt',chr)
	assign(obj, dt)
	dtlist[[j]] <- get(obj)
	j <- j + 1
}

master <- rbindlist(dtlist)
write.table(master, '/mnt/lab_data/montgomery/nicolerg/rfmix-outputs/rfmix.viterbi.master.tsv', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

master[1:10,1:10]
