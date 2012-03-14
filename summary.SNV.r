stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==22)
	{
	gene=as.character(stdin[6])
	input.nonsense=as.character(stdin[7])
	input.missense=as.character(stdin[8])
	input.codingsynonymous=as.character(stdin[9])
	input.codingnotMod3=as.character(stdin[10])
	input.splice3=as.character(stdin[11])
	input.splice5=as.character(stdin[12])
	input.utr3=as.character(stdin[13])
	input.utr5=as.character(stdin[14])
	output.nonsense=as.character(stdin[15])
	output.missense=as.character(stdin[16])
	output.codingsynonymous=as.character(stdin[17])
	output.codingnotMod3=as.character(stdin[18])
	output.splice3=as.character(stdin[19])
	output.splice5=as.character(stdin[20])
	output.utr3=as.character(stdin[21])
	output.utr5=as.character(stdin[22])
	}
}
gene=as.matrix(read.table(file=gene))
nonsense=as.matrix(read.table(file=input.nonsense))
missense=as.matrix(read.table(file=input.missense))
codingsynonymous=as.matrix(read.table(file=input.codingsynonymous))
codingnotMod3=as.matrix(read.table(file=input.codingnotMod3))
splice3=as.matrix(read.table(file=input.splice3))
splice5=as.matrix(read.table(file=input.splice5))
utr3=as.matrix(read.table(file=input.utr3))
utr5=as.matrix(read.table(file=input.utr5))
out.nonsense=matrix(0,nrow=nrow(gene), ncol=2)
out.missense=matrix(0,nrow=nrow(gene), ncol=2)
out.codingsynonymous=matrix(0,nrow=nrow(gene), ncol=2)
out.codingnotMod3=matrix(0,nrow=nrow(gene), ncol=2)
out.splice3=matrix(0,nrow=nrow(gene), ncol=2)
out.splice5=matrix(0,nrow=nrow(gene), ncol=2)
out.utr3=matrix(0,nrow=nrow(gene), ncol=2)
out.utr5=matrix(0,nrow=nrow(gene), ncol=2)

out.nonsense[,1] = out.missense[,1] = out.codingsynonymous[,1] = out.codingnotMod3[,1] = out.splice3[,1] = out.splice5[,1] = out.utr3[,1] = out.utr5[,1] =as.character(gene[,1])

for (i in 1:nrow(nonsense))
{
	out.nonsense[out.nonsense[,1]==nonsense[i],2] = as.numeric(out.nonsense[out.nonsense[,1]==nonsense[i],2]) + 1
}
for (i in 1:nrow(missense))
{
	out.missense[out.missense[,1]==missense[i],2] = as.numeric(out.missense[out.missense[,1]==missense[i],2]) + 1
}
for (i in 1:nrow(codingsynonymous))
{
	out.codingsynonymous[out.codingsynonymous[,1]==codingsynonymous[i],2] = as.numeric(out.codingsynonymous[out.codingsynonymous[,1]==codingsynonymous[i],2]) + 1
}
for (i in 1:nrow(codingnotMod3))
{
	out.codingnotMod3[out.codingnotMod3[,1]==codingnotMod3[i],2] = as.numeric(out.codingnotMod3[out.codingnotMod3[,1]==codingnotMod3[i],2]) + 1
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

write.table(out.nonsense, file = output.nonsense, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.missense, file = output.missense, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.codingsynonymous, file = output.codingsynonymous, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.codingnotMod3, file = output.codingnotMod3, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.splice3, file = output.splice3, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.splice5, file = output.splice5, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.utr3, file = output.utr3, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.utr5, file = output.utr5, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
