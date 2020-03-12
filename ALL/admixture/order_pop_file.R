library(data.table)

args <- commandArgs(trailingOnly=T)
fam_file <- args[1]
tg <- args[2]

fam <- fread(fam_file,sep=' ',header=F)
tg <- fread(tg,sep='\t',header=T)

fam <- fam[,.(V1)]
m <- merge(fam, tg, by.x='V1',by.y='Subject',all.x=T)
m[is.na(m)] <- '-'

head(m)
m <- m[match(fam[,V1], m[,V1])]
head(m)
nrow(m)
nrow(fam)

outfile <- gsub('fam','pop.txt',fam_file)
write.table(m, outfile, col.names=F,row.names=F,quote=F,sep='\t')
