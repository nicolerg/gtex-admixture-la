df <- read.delim("~/Documents/Surya/gtex-admix/plots/geuvadis/peer-eqtl-counts-fixed.tsv", sep='\t', header=TRUE)
library(ggplot2)
pdf("~/Documents/Surya/gtex-admix/plots/geuvadis/peer-eqtl-discovery-chr22.pdf",width=8, height=6)
ggplot(df, aes(x=K))+
  geom_line(aes(y=N_eur, colour="N_eur")) +
  geom_line(aes(y=N_afr, colour="N_afr")) +
  geom_vline(xintercept=13, linetype="dashed") +
  scale_colour_manual(labels=c("N_eur"="EUR","N_afr"="AFR"), values=c(N_eur="blue",N_afr="#FFB933"), name="Population") +
  labs(x="N PEER factors",y="N cis-eQTL genes on chr22") 
dev.off()
