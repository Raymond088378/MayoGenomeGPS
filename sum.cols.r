stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==7)
	{
	in.summary=as.character(stdin[6])
	output.sum=as.character(stdin[7])
	}
}
in.summary=read.table(file=in.summary)
out.sum=matrix(0,nrow=nrow(in.summary), ncol=1)

for (i in 1:nrow(in.summary))
{
	out.sum[i] = rowSums(in.summary[i,])
}

write.table(out.sum, file = output.sum, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
