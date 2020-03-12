# make a bar plot of the fraction of f2 variants with matching GTEx LA and 1000 Genomes populations for each chr

library(ggplot2)
library(reshape)

## in bash:
# sed 's/.* //' 1kg-af-ancestry-crossref.out | sed '/^0/!d' | pr -ts" " --columns 22 > f2-1kg-frac-overlap-bychr.tsv
# echo "FRAC" >> tmp; echo "AFR" >> tmp; echo "ASN" >> tmp; echo "EUR" >> tmp
# echo "chr1	chr2	chr3	chr4	chr5	chr6	chr7	chr8	chr9	chr10	chr11	chr12	chr13	chr14	chr15	chr16	chr17	chr18	chr19	chr20	chr12	chr22" >> tmp2
# cat f2-1kg-frac-overlap-bychr.tsv >> tmp2 
# rm f2-1kg-frac-overlap-bychr.tsv
# paste tmp tmp2 > f2-1kg-frac-overlap-bychr.tsv

f <- '/mnt/lab_data/montgomery/nicolerg/f2-variants/4_1kgAF-gtexLA-match/f2-1kg-frac-overlap-bychr.tsv'

frac <- read.table(f, sep='\t', header=TRUE)
melted <- melt(frac, id='FRAC')
head(melted)

pdf(file='/users/nicolerg/gtex-admix/f2-1kg-frac-overlap.pdf',width=12,height=5)
ggplot(melted, aes(fill=FRAC, y=value, x=variable)) + 
    geom_bar(position="dodge", stat="identity", colour="black",linetype="solid") +
    labs(y="Matching Fraction",x="Chromosome") +
    scale_fill_manual(name="Population", values=c(ASN="#ADFF33",EUR="#0000FF",AFR="#FF9900",TOTAL="gray"))
dev.off()
