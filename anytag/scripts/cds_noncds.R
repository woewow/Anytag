df<-read.table("all.cds.noncds.stat2")
pdf("cds_stat.pdf", width=8)
barplot(t(as.matrix(df)), beside=T, names.arg=rownames(df), cex.names=0.7, legend=names(df), xlab="Coverage", ylab="Frequency", col = colors()[c(36,448,275,91)], main="CS reads and sanger reads in CDS and nonCDS distribution", args.legend=list(x="topright"))
dev.off()
