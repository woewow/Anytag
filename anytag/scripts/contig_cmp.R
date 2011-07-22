long <- read.table("long.contig.length")
short <- read.table("short.contig.length")
sanger <- read.table("sanger.contig.length")
#long1 <- log10(long$V1)
#short1 <- log10(short$V1)
#sanger1 <- log10(sanger$V1)
#pdf("contig_cmp.pdf")
par(mfrow=c(3,1))
long1 <- seq(0, max(long$V1), by=1000)
short1 <- seq(0, max(short$V1), by=1000)
sanger1 <- seq(0, max(sanger$V1), by=1000)
j <- 0
for(i in seq(0, max(long$V1), by=1000)) {
	long1[j] <- sum(long$V1[long$V1>=i&long$V1<i+1000])
	j <- j+1
}
j <- 0
for(i in seq(0, max(short$V1), by=1000)) {
	short1[j] <- sum(short$V1[short$V1>=i&short$V1<i+1000])
	j <- j+1
}
j <- 0
for(i in seq(0, max(sanger$V1), by=1000)) {
	sanger1[j] <- sum(sanger$V1[sanger$V1>=i&sanger$V1<i+1000])
	j <- j+1
}
plot(long1/sum(long1), type="l", ylab="", las=2, ylim=c(0,.15),xlim=c(0,160),xaxt="n",xlab="")
abline(v=18119/1000, col="green",lwd=2)
abline(v=2824/1000, col="blue",lwd=2)
abline(v=6152.6/1000, col="red",lwd=2)
title(xlab="PCAP assembled contig length of CS reads", ylab="Fraction of bases", cex.lab=1.2)
legend("topright", c("N50","N90", "Mean"), lwd=c(2,2,2),lty=c(1,1, 1), col=c("green","blue", "red"))
axis(1, las=1, labels=seq(0,300, by=10)*1000, at=seq(0,300,by=10))
plot(short1/sum(short1), type="l", ylab="", las=2,ylim=c(0,.15),xlim=c(0,160),xaxt="n", xlab="")
abline(v=6389/1000, col="green",lwd=2)
abline(v=958/1000, col="blue",lwd=2)
abline(v=1963/1000, col="red",lwd=2)
title(xlab="Velvet assembled contig length of short reads", ylab="Fraction of bases", cex.lab=1.2)
legend("topright", c("N50","N90", "Mean"), lwd=c(2,2,2),lty=c(1,1, 1), col=c("green","blue", "red"))
axis(1, las=1, labels=seq(0,300, by=10)*1000, at=seq(0,300,by=10))
plot(sanger1/sum(sanger1), type="l", ylab="", las=2,ylim=c(0,.15),xlim=c(0,160),xaxt="n", xlab="")
abline(v=29351/1000, col="green",lwd=2)
abline(v=1014/1000, col="blue",lwd=2)
abline(v=4523.9/1000, col="red",lwd=2)
title(xlab="PCAP assembled contig length of sanger reads", ylab="Fraction of bases", cex.lab=1.2)
legend("topright", c("N50","N90", "Mean"), lwd=c(2,2,2),lty=c(1,1, 1), col=c("green","blue", "red"))
axis(1, las=1, labels=seq(0,300, by=10)*1000, at=seq(0,300,by=10))

#hist((long$V1), breaks=seq(0,max(long$V1)+1000,by=1000),  xlim=c(0, 150000), ylim=c(0, .002))
#hist((short$V1), breaks=seq(0,max(short$V1)+1000,by=1000),  xlim=c(0, 150000), ylim=c(0, .002))
#hist((sanger$V1), breaks=seq(0,max(sanger$V1)+1000,by=1000),  xlim=c(0, 150000), ylim=c(0, .002))

#hist(log2(long$V1), breaks=seq(0,log2(max(long$V1))+log2(100),by=log2(100)), freq=F, ylim=c(0, .6))
#hist(log2(short$V1), breaks=seq(0,log2(max(short$V1))+log2(100),by=log2(100)), freq=F, ylim=c(0, .2))
#hist(log2(sanger$V1), breaks=seq(0,log2(max(sanger$V1))+log2(100),by=log2(100)), freq=F, ylim=c(0, .2))

#hist(long1, breaks=c(0, log10(seq(1000,max(long$V1)+1000,by=1000))), freq=T, xlim=c(0, 20), xlab="PCAP assembled contig length (log2) of CS reads", main="", ylim=c(0,12000), lwd=4)
#hist(short1, breaks=c(0, log10(seq(1000,max(short$V1)+1000,by=1000))),freq=T,  xlim=c(0,20), xlab="Velvet assembled contig length (log2) of short reads", main="", ylim=c(0,12000),lwd=4)
#hist(sanger1,breaks=c(0, log10(seq(1000,max(sanger$V1)+100,by=1000))),freq=T,xlim=c(0,20), xlab="PCAP assembled contig length (log2) of sanger reads", main="", ylim=c(0,12000),lwd=4)

#dev.off()