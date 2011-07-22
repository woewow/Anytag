df<-read.table("all.repeat.nonrep.stat2")
pdf("repeat_stat.pdf", width=10)
barplot(t(as.matrix(df)), beside=T, names.arg=rownames(df), cex.names=0.7, legend=names(df), xlab="Coverage", ylab="Frequency", col = colors()[c(113,140,84,490)], main="CS reads and sanger reads in repeat and nonrepeat distribution", args.legend=list(x="topright"))
dev.off()
