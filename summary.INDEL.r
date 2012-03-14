stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==18)
	{
	gene=as.character(stdin[6])
	input.coding=as.character(stdin[7])
	input.frameshift=as.character(stdin[8])
	input.splice3=as.character(stdin[9])
	input.splice5=as.character(stdin[10])
	input.utr3=as.character(stdin[11])
	input.utr5=as.character(stdin[12])
	output.coding=as.character(stdin[13])
	output.frameshift=as.character(stdin[14])
	output.splice3=as.character(stdin[15])
	output.splice5=as.character(stdin[16])
	output.utr3=as.character(stdin[17])
	output.utr5=as.character(stdin[18])
	}
}
gene=as.matrix(read.table(file=gene))
coding=as.matrix(read.table(file=input.coding))
frameshift=as.matrix(read.table(file=input.frameshift))
splice3=as.matrix(read.table(file=input.splice3))
splice5=as.matrix(read.table(file=input.splice5))
utr3=as.matrix(read.table(file=input.utr3))
utr5=as.matrix(read.table(file=input.utr5))
out.coding=matrix(0,nrow=nrow(gene), ncol=2)
out.frameshift=matrix(0,nrow=nrow(gene), ncol=2)
out.splice3=matrix(0,nrow=nrow(gene), ncol=2)
out.splice5=matrix(0,nrow=nrow(gene), ncol=2)
out.utr3=matrix(0,nrow=nrow(gene), ncol=2)
out.utr5=matrix(0,nrow=nrow(gene), ncol=2)

out.coding[,1] = out.frameshift[,1] = out.splice3[,1] = out.splice5[,1] = out.utr3[,1] = out.utr5[,1] =as.character(gene[,1])

for (i in 1:nrow(coding))
{
	out.coding[out.coding[,1]==coding[i],2] = as.numeric(out.coding[out.coding[,1]==coding[i],2]) + 1
}
for (i in 1:nrow(frameshift))
{
	out.frameshift[out.frameshift[,1]==frameshift[i],2] = as.numeric(out.frameshift[out.frameshift[,1]==frameshift[i],2]) + 1
}
for (i in 1:nrow(splice3))
{
	out.splice3[out.splice3[,1]==splice3[i],2] = as.numeric(out.splice3[out.splice3[,1]==splice3[i],2]) + 1
}
for (i in 1:nrow(splice5))
{
	out.splice5[out.splice5[,1]==splice5[i],2] = as.numeric(out.splice5[out.splice5[,1]==splice5[i],2]) + 1
}
for (i in 1:nrow(utr3))
{
	out.utr3[out.utr3[,1]==utr3[i],2] = as.numeric(out.utr3[out.utr3[,1]==utr3[i],2]) + 1
}
for (i in 1:nrow(utr5))
{
	out.utr5[out.utr5[,1]==utr5[i],2] = as.numeric(out.utr5[out.utr5[,1]==utr5[i],2]) + 1
}

write.table(out.coding, file = output.coding, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.frameshift, file = output.frameshift, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.splice3, file = output.splice3, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.splice5, file = output.splice5, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.utr3, file = output.utr3, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.utr5, file = output.utr5, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
