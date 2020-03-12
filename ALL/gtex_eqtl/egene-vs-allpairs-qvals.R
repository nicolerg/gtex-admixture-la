df <- read.delim("/users/nicolerg/gtex-admix/metadata/intersect_qval.tsv",sep='\t',header=FALSE)
colnames(df) <- c("gene_id","stat","q_val")

df$q_val <- -log(as.numeric(df$q_val))
df_match <- df[as.character(df$stat)=="match",]
df_uniq <- df[as.character(df$stat)=="none",]

jpeg("/users/nicolerg/gtex-admix/qval-unique.jpg", width=1600, height=1400, res=200)
hist(as.numeric(df_uniq$q_val),breaks=50,main=paste("Q Values for Unique Genes (N = ",nrow(df_uniq),")",sep=''),xlab="Q Value (-log)",ylim=c(0,600))
dev.off()

jpeg("/users/nicolerg/gtex-admix/qval-match.jpg", width=1600, height=1400, res=200)
hist(as.numeric(df_match$q_val),breaks=50,main=paste("Q Values for Matching Genes (N = ",nrow(df_match),")",sep=''),xlab="Q Value (-log)",ylim=c(0,4000))
dev.off()