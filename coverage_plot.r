#; Input file name should be of the form sample1.txt, sample2.txt...sample8.txt. The files should contain the %values at each depth of coverage.
#; Usage: Rscript coverage_plot.r [target region] [sample names space seperated]
#;Example: Rscript coverage_plot.r 50377064 A244 A345

stdin <- commandArgs(TRUE) 

linecolors <- c("blue","green","red","black","brown","violet","gray50","orange","olivedrab","cyan","royalblue","yellow","darkorchid4","wheat4","seagreen","gray25","darkslateblue","darkgoldenrod4","cadetblue","chartreuse3","violetred4","palegreen4","skyblue4","pink4")
symbol = 1:50
if(length(stdin) > 1){
	target <- as.integer(stdin[1])
	samplenames <- c()
	jpeg(file="coverage.jpeg", width=1190, height=1190, quality=100, res=100)
	for (i in 1:(length(stdin)-1)){
		if(i <= length(linecolors)){
			samplename <- as.character(stdin[i+1])
			samplenames <- cbind(samplenames,samplename)
			coveragefile <- paste(samplename,".coverage.out", sep="", collapse=NULL)
			samplematrix <- as.matrix(read.table(file=coveragefile))
			percent <- (samplematrix/target)*100
			linecolor <- linecolors[i]
			sym <- i
			if(i == 1){
				plot(percent,xlim=c(1,nrow(samplematrix)),ylim=c(0,115),main="Coverage across target regions at different depths of coverage",xlab="Depth of coverage",ylab="Percent coverage", pch=sym,col=linecolor,type="o",cex=0.6,lwd=2)
			}else{
				lines(percent,col=linecolor,pch=sym,type="o",cex=0.6,lwd=2)
			}
		}
	}
	legend(x="bottomleft",samplenames,inset=0.02,cex=0.7,col=linecolors[1:length(samplenames)],pch=symbol[1:length(samplenames)],lty=1,lwd=2)
	dev.off()
	
}else{
	print("Did not supply enough arguments to plot coverage!")
}
