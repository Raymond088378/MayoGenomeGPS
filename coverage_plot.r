#; Input file name should be of the form sample1.txt, sample2.txt...sample8.txt. The files should contain the %values at each depth of coverage.
#; Usage: Rscript coverage_plot.r [target region] [sample names space seperated]
#;Example: Rscript coverage_plot.r 50377064 A244 A345
stdin <- commandArgs(TRUE) 

if(length(stdin) > 1){
	jpeg(file="coverage.jpeg", width=1190, height=1190, quality=100, res=100)
	target <- as.integer(stdin[1])
	sp.max <- length(stdin)-1
	if (sp.max >1 )	{
		palette(rainbow(sp.max))
	}
	all <- c()
	samplenames <- c()
	for (i in 1:sp.max) {
		samplename <- as.character(stdin[i+1])
		samplenames <- cbind(samplenames,samplename)
		coveragefile <- paste(samplename,".coverage.out", sep="", collapse=NULL)
		samplematrix <- as.matrix(read.table(file=coveragefile))
		percent <- (samplematrix/target)*100
		all <- cbind(all, percent)
		linecolor <- i
		sym <- i
		if(i == 1){
			plot(percent,xlim=c(1,nrow(samplematrix)),ylim=c(0,quantile(percent,0.95)),main="Coverage across target regions at different depths of coverage",xlab="Depth of coverage",ylab="Percent coverage", pch=sym,col=linecolor,type="o",cex=0.6,lwd=2)
		}else{
			lines(percent,col=linecolor,pch=sym,type="o",cex=0.6,lwd=2)
		}
	}
	mean.percent <- rowMeans(all)
	mean.len <- length(mean.percent)
	sd1 <- c()
	sd2 <- c()
	for (i in 1:mean.len) {
		row.tmp <- all[i,]
		row.sd <- sd(row.tmp)
		row.sd1 <- mean.percent[i] - 3*row.sd	
		row.sd2 <- mean.percent[i] + 3*row.sd	
		sd1 <- c(sd1, row.sd1) 
		sd2 <- c(sd2, row.sd2)
	}
	linecolors <- seq(1:sp.max)
	symbol <- seq(1:sp.max) 
	
	legend(x="bottomleft",samplenames,inset=0.02,cex=0.7,col=linecolors,pch=symbol,lty=1,lwd=2)
	dev.off()
	
}else{
	print("Did not supply enough arguments to plot coverage!")
}
