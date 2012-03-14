stdin = commandArgs() 
## Input file 1: values for used reads, mapped reads, genome mapped and junction mapped
## Input file 2: samples names
for (i in 1:1)
{
	if (length(stdin)==7)
	{
	values=as.character(stdin[6])
	samplename=as.character(stdin[7])
	}
}
data=read.table(file=values, header=T, sep="\t")
samples=as.matrix(read.table(file=samplename))
max_y=max(data)
min_y=min(data)
num=nrow(samples)
plot_colors <- c("blue","red","forestgreen", "black")
#pdf(file="QC_ReadsDistribution.pdf", width=20, heigh=10)
png(filename="QC_ReadsDistribution.png", bg="white", res=100)
plot(data$UsedReads, type="o", col=plot_colors[1], lwd=1, ylim=c(min_y,max_y), axes=FALSE, ann=FALSE)
#axis(1, at=1:num, lab=F)
#text(axTicks(1), par("usr")[3] - 2, srt=45, adj=1,labels=c(samples),xpd=T, cex=0.8)
#text(8*(1:num), srt=45, par("usr")[3] - 2, adj=1, labels=c(samples), xpd=T, cex=0.4)
text(1:num, srt=45, par("usr")[3] - 2, adj=1, labels=c(samples), xpd=T, cex=0.7)

axis(2, las=1, cex.axis=0.8)
box()
lines(data$MappedReads, type="o", pch=22, lty=2, lwd=1, col=plot_colors[2], ylim=c(min_y,max_y))
lines(data$GenomeMapped, type="o", pch=23, lty=3, lwd=1, col=plot_colors[3], ylim=c(min_y,max_y))
lines(data$JunctionMapped, type="o", pch=24, lty=4, lwd=1, col=plot_colors[4], ylim=c(min_y,max_y))
title(main="QC: Distribution of Reads", col.main="red", font.main=4)
title(xlab= "Samples", col.lab=rgb(0,0.5,0))
title(ylab= "Counts", col.lab=rgb(0,0.5,0))
par(xpd=TRUE)
# if (num > 8)
# {
	# legend(450, (max_y/3), names(data), cex=0.8, col=plot_colors, lty=1:4, lwd=1, bty="n")
# }else{
	# legend((num-2), (max_y/2), names(data), cex=0.5, col=plot_colors, lty=1:4, lwd=1, bty="n")
# }
legend((num-2), (max_y/2), names(data), cex=0.5, col=plot_colors, lty=1:4, lwd=1, bty="n")
dev.off()


