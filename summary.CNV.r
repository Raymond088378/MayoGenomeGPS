stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==10)
	{
	gene=as.character(stdin[6])
	input.DEL=as.character(stdin[7])
	input.DUP=as.character(stdin[8])
	output.DEL=as.character(stdin[9])
	output.DUP=as.character(stdin[10])
	}
}
gene=as.matrix(read.table(file=gene))
DEL=as.matrix(read.table(file=input.DEL))
DUP=as.matrix(read.table(file=input.DUP))
out.DEL=matrix(0,nrow=nrow(gene), ncol=2)
out.DUP=matrix(0,nrow=nrow(gene), ncol=2)

out.DEL[,1] = out.DUP[,1] = as.character(gene[,1])

for (i in 1:nrow(DEL))
{
	out.DEL[out.DEL[,1]==DEL[i,1],2] = as.numeric(out.DEL[out.DEL[,1]==DEL[i,1],2]) + 1
}
for (i in 1:nrow(DUP))
{
	out.DUP[out.DUP[,1]==DUP[i,1],2] = as.numeric(out.DUP[out.DUP[,1]==DUP[i,1],2]) + 1
}

write.table(out.DEL, file = output.DEL, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.DUP, file = output.DUP, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
