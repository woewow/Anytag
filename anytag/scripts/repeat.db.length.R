df<-read.table("repeat.db.length")
pdf("repeat.db.length.distribution.pdf")
hist(df$V1[df$V1<1000], breaks=100, main="repeat db length(<1000bp) distribution", xlab="length")
dev.off()
