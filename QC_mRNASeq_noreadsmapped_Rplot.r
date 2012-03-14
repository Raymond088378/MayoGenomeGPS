stdin = commandArgs() 
## Input file 1: number of genes with no counts
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
plot_colors <- c("blue")
#pdf(file="QC_ReadsDistribution.pdf", width=20, heigh=10)
png(filename="QC_AllSamples_GeneCount_NoReadsMapped.png", bg="white", res=100)
plot(data$Count, type="o", col=plot_colors[1], ylim=c(min_y,max_y), lwd=1, axes=FALSE, ann=FALSE)
text(1:num, srt=45, par("usr")[3] - 2, adj=1, labels=c(samples), xpd=T, cex=0.7)
axis(2, las=1, cex.axis=0.8)
box()
title(main="QC: Number of genes with no reads mapped", col.main="red", font.main=4)
title(xlab= "Samples", col.lab=rgb(0,0.5,0))
title(ylab= "Counts", col.lab=rgb(0,0.5,0))
par(xpd=TRUE)
dev.off()


