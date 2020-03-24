library(data.table)

args <- commandArgs(trailing=T)
superpops <- unname(unlist(strsplit(args[1],',')))
genotyped <- fread(args[2],sep='\t',header=F) # samples in merged VCFs
tg <- fread(args[3],sep='\t',header=T) # 1000 Genomes population codes
gtex <- fread(args[4],sep='\t',header=F) # list of GTEx v8 admixed individuals (gtex-admixed0.9.txt)

head(genotyped)

tg <- tg[Subject %in% genotyped[,V1]]

# remove admixed subpopulations
h_ad <- c('ASW','ACB')
tg <- tg[! Subpop %in% h_ad]

head(tg)

# # select 100 reference panels from each superpopulation (equal numbers from each subpopulation)
# keep <- c()
# for (pop in superpops){
# 	n_subpops = length(unique(tg[Superpop==pop,Subpop]))
# 	for (sub in unique(tg[Superpop==pop,Subpop])){
# 		selected <- sample(tg[Superpop==pop&Subpop==sub,Subject], round(100/n_subpops), replace=F)
# 		keep <- c(keep, selected)
# 	}
# }

# select all reference panels from each superpopulation
keep = tg[Superpop %in% superpops,Subject]

# prepend GTEx IDs
keep <- c(gtex[,V1], keep)

# keep populations handy
keep <- data.table(SUBJID=keep)
keep <- merge(keep, tg, by.x='SUBJID', by.y='Subject', all.x=T)
keep[is.na(keep)] <- '-'

write.table(keep, 'keep_subjects.txt', sep='\t', col.names=F, row.names=F, quote=F)
