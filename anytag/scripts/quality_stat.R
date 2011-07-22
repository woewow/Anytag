#df <- read.table("shortread.qual.count.txt2")
#df2 <- read.table("longread.qual.count.txt")
df <- read.table("sample.qual")
dferr <- read.table("sample.err.qual.count")
df2 <- read.table("sample.qual.count.longread")
df2err <- read.table("sample.err.qual.count.longread")
df3 <- read.table("sample.qual.count.sanger")
df3err <- read.table("sample.err.qual.count.sanger")

freq <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	freq[j] <- sum(df$V1[(df$V3<i)&(df$V3>=i-10)])
	j <- j+1 
}
err <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	err[j] <- sum(dferr$V1[dferr$V2<i&dferr$V2>=i-10])
	j <- j+1
}
cnt <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	cnt[j] <- sum(df$V2[df$V3<i&df$V3>=i-10])
	j <- j+1
}

#print(freq)
freq2 <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	freq2[j] <- sum(df2$V1[df2$V2<i&df2$V2>=i-10])
	j <- j+1
}
err2 <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	err2[j] <- sum(df2err$V1[df2err$V2<i&df2err$V2>=i-10])
	j <- j+1
}
cnt2 <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	cnt2[j] <- sum(df2$V3[df2$V2<i&df2$V2>=i-10])
	j <- j+1
}

freq3 <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	freq3[j] <- sum(df3$V1[df3$V2<i&df3$V2>=i-10])
	j <- j+1
}
err3 <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	err3[j] <- sum(df3err$V1[df3err$V2<i&df3err$V2>=i-10])
	j <- j+1
}
cnt3 <- seq_len(10)
j <- 1
for (i in seq(10,100,by=10)) {
	cnt3[j] <- sum(df3$V3[df3$V2<i&df3$V2>=i-10])
	j <- j+1
}

#print(freq2)
pdf("quality_eval.pdf",width=10)
par(mar=c(5, 4, 4, 6) + 0.1)
plot(freq*100,type="b", lty=3, lwd=2, pch=20, ylim=c(0,100), col=colors()[76],ylab="Percentage", xlab="Quality(10 based intervals)",axes=F)
points(freq2*100, type="b", lwd=2, lty=3, pch=4,col=colors()[26])
points(freq3*100, type="b", lwd=2, lty=3, pch=2,col=colors()[51])
axis(1,labels=c(paste(seq(0,80,by=10),9,sep="-"),">=90"),at=seq(1,10))
axis(2,las=1)
box()
par(new = T)
plot(log(err2/cnt2), type="b",ylim=c(-1,-10),lwd=2, lty=1, pch=4, col=colors()[26], xlab="",ylab="", axes=F)
short_err <- err/cnt
for(j in which(is.na(short_err))) {
	short_err[j] <- 0
}
short_err3 <- err3/cnt3
for(j in which(is.na(short_err))) {
	short_err3[j] <- 0
}
points(log(short_err), type="b", ylim=c(-1,-10), lwd=2, lty=1, pch=20, col=colors()[76])
points(log(short_err3), type="b", ylim=c(-1,-10), lwd=2, lty=1, pch=2, col=colors()[51])
mtext("Error rate",side=4, line=3.5)
axis(4,las=1,labels=c(paste("10e",-1:-10,sep="")),at=(-1:-10))
legend("topleft", legend=c("CS reads percent bases", "CS reads error rate", "short reads percent bases", "short reads error rate", "sanger reads percent bases", "sanger reads error rate"), col=colors()[c(26,26,76,76,51,51)], pch=c(20,20,4,4,2,2), lty=c(3,1,3,1,3,1), lwd=c(2,2,2,2,2,2),text.col=c("gray","gray","gray","gray","gray","gray"))
dev.off()

