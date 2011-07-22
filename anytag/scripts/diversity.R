df <- read.table("longread.heterozygosity")
pdf("diversity_in_genomic_region.pdf")
#tiff(filename = "diversity_in_genomic_region.tif",width = 960, height = 960, res=1000)
#dp <- barplot(df$V1*100, names.arg=df$V3,las=3,ylab="Diversity (%)", col="dark blue",main="diversity between w1118 and dm3 in different genomic regions")
dp <- barplot(df$V1*100, axes=F,axisnames=F, ylab="Diversity (%)",main="diversity between w1118 and dm3 in different genomic regions", col="purple2")
text(dp, par("usr")[3], labels=df$V3, srt=45, adj=c(1.1,1.1), xpd=T, cex=.9)
axis(2, las = 1)
dev.off()

