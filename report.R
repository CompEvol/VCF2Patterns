
library("ComMA")
setwd("~/WorkSpace/VCF2Patterns/result")

df.mu <- ComMA::readFile("snp_frequencies.txt", header=FALSE, row.names=NULL) 
colnames(df.mu) <- c("Counts", "Samples")
# folded site frequency spectrum
df.mu[df.mu[,1]>46,1] <- 92 - df.mu[df.mu[,1]>46,1]
df.mu.aggr <- aggregate(Samples ~ Counts, data = df.mu, sum)
bar.chart <- ggBarChart(df.mu.aggr, x.id="Counts", y.id="Samples", y.scale="log", title="SNPs Count Across Samples", y.lim.cart = c(10000, NA))
pdfGgplot(bar.chart, fig.path="snps.pdf")

df.idel <- ComMA::readFile("indel_frequencies.txt", header=FALSE, row.names=NULL) 
colnames(df.idel) <- c("Counts", "Samples")
# folded site frequency spectrum
df.idel[df.idel[,1]>46,1] <- 92 - df.idel[df.idel[,1]>46,1]
df.idel.aggr <- aggregate(Samples ~ Counts, data = df.idel, sum)
bar.chart <- ggBarChart(df.idel.aggr, x.id="Counts", y.id="Samples", y.scale="log", title="Indels Count Across Samples", y.lim.cart = c(100, NA))
pdfGgplot(bar.chart, fig.path="idels.pdf")

setwd("~/WorkSpace/VCF2Patterns/result/modern-ss")

df.mu <- ComMA::readFile("snp_frequencies.txt", header=FALSE, row.names=NULL) 
colnames(df.mu) <- c("Counts", "Samples")
# folded site frequency spectrum
df.mu[df.mu[,1]>46,1] <- 92 - df.mu[df.mu[,1]>46,1]
df.mu.aggr <- aggregate(Samples ~ Counts, data = df.mu, sum)
bar.chart <- ggBarChart(df.mu.aggr, x.id="Counts", y.id="Samples", y.scale="log", title="SNPs Count Across Samples (Modern)", y.lim.cart = c(10000, NA))
pdfGgplot(bar.chart, fig.path="snps-modern.pdf")

df.idel <- ComMA::readFile("indel_frequencies.txt", header=FALSE, row.names=NULL) 
colnames(df.idel) <- c("Counts", "Samples")
# folded site frequency spectrum
df.idel[df.idel[,1]>46,1] <- 92 - df.idel[df.idel[,1]>46,1]
df.idel.aggr <- aggregate(Samples ~ Counts, data = df.idel, sum)
bar.chart <- ggBarChart(df.idel.aggr, x.id="Counts", y.id="Samples", y.scale="log", title="Indels Count Across Samples (Modern)", y.lim.cart = c(100, NA))
pdfGgplot(bar.chart, fig.path="idels-modern.pdf")


