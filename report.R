
library("ComMA")
setwd("~/WorkSpace/VCF2Patterns/result")

df.mu <- ComMA::readFile("snp_frequencies.txt", header=FALSE, row.names=NULL) 
colnames(df.mu) <- c("Counts", "Samples")
df.mu[df.mu[,1]>46,1] <- 92 - df.mu[df.mu[,1]>46,1]
df.mu.aggr <- aggregate(Samples ~ Counts, data = df.mu, sum)
bar.chart <- ggBarChart(df.mu.aggr, x.id="Counts", y.id="Samples", y.scale="log", title="SNPs Count Across Samples")
pdfGgplot(bar.chart, fig.path="snps.pdf")

df.idel <- ComMA::readFile("indel_frequencies.txt", header=FALSE, row.names=NULL) 
colnames(df.idel) <- c("Counts", "Samples")
df.idel[df.idel[,1]>46,1] <- 92 - df.idel[df.idel[,1]>46,1]
df.idel.aggr <- aggregate(Samples ~ Counts, data = df.idel, sum)
bar.chart <- ggBarChart(df.idel.aggr, x.id="Counts", y.id="Samples", y.scale="log", title="Indels Count Across Samples")
pdfGgplot(bar.chart, fig.path="idels.pdf")
