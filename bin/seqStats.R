#!/usr/bin/Rscript

#R CMD BATCH ./bin/seqStats.R

seqStats <- read.table("data/sequences/summary-stats.tsv", sep="\t", header=TRUE)
pdf(file="./docs/figures/seqStats.pdf", width=18, height=18)
par(cex=2.0, las=2,mar=c(4, 5, 4, 2) + 0.1)
plot(seqStats$CG.content, log10(seqStats$Length), xlim=c(0.1,0.9), ylim=c(1,9), xlab="C+G content", ylab="Sequence length",main="",yaxt = "n", pch=21, col="black", bg="magenta")
#grid
xp <- seq(0.1,0.9, by=0.1)
yp <- 1:9
for(i in 1:length(xp)){
	lines(c(xp[i],xp[i]), c(0,10))
}
for(j in 1:length(yp)){
	lines(c(0,1), c(yp[j],yp[j]))
}
tcks<-10^(1:9); axis(2,at=log10(tcks), tcks)
dev.off()


