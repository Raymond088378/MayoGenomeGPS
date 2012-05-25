stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==32)
	{
	input.gene=as.character(stdin[6])
	input.EXON_DELETED=as.character(stdin[7])
	input.FRAME_SHIFT=as.character(stdin[8])
	input.CODON_CHANGE=as.character(stdin[9])
	input.UTR_5_DELETED=as.character(stdin[10])
	input.UTR_3_DELETED=as.character(stdin[11])
	input.CODON_INSERTION=as.character(stdin[12])
	input.CODON_CHANGE_PLUS_CODON_INSERTION=as.character(stdin[13])
	input.CODON_DELETION=as.character(stdin[14])
	input.CODON_CHANGE_PLUS_CODON_DELETION=as.character(stdin[15])
	input.SPLICE_SITE_ACCEPTOR=as.character(stdin[16])
	input.SPLICE_SITE_DONOR=as.character(stdin[17])
	input.UTR_5_PRIME=as.character(stdin[18])
	input.UTR_3_PRIME=as.character(stdin[19])
	output.EXON_DELETED=as.character(stdin[20])
	output.FRAME_SHIFT=as.character(stdin[21])
	output.CODON_CHANGE=as.character(stdin[22])
	output.UTR_5_DELETED=as.character(stdin[23])
	output.UTR_3_DELETED=as.character(stdin[24])
	output.CODON_INSERTION=as.character(stdin[25])
	output.CODON_CHANGE_PLUS_CODON_INSERTION=as.character(stdin[26])
	output.CODON_DELETION=as.character(stdin[27])
	output.CODON_CHANGE_PLUS_CODON_DELETION=as.character(stdin[28])
	output.SPLICE_SITE_ACCEPTOR=as.character(stdin[29])
	output.SPLICE_SITE_DONOR=as.character(stdin[30])
	output.UTR_5_PRIME=as.character(stdin[31])
	output.UTR_3_PRIME=as.character(stdin[32])

	}
}
gene=as.matrix(read.table(file=input.gene))
EXON_DELETED=as.matrix(read.table(file=input.EXON_DELETED))
FRAME_SHIFT=as.matrix(read.table(file=input.FRAME_SHIFT))
CODON_CHANGE=as.matrix(read.table(file=input.CODON_CHANGE))
UTR_5_DELETED=as.matrix(read.table(file=input.UTR_5_DELETED))
UTR_3_DELETED=as.matrix(read.table(file=input.UTR_3_DELETED))
CODON_INSERTION=as.matrix(read.table(file=input.CODON_INSERTION))
CODON_CHANGE_PLUS_CODON_INSERTION=as.matrix(read.table(file=input.CODON_CHANGE_PLUS_CODON_INSERTION))
CODON_DELETION=as.matrix(read.table(file=input.CODON_DELETION))
CODON_CHANGE_PLUS_CODON_DELETION=as.matrix(read.table(file=input.CODON_CHANGE_PLUS_CODON_DELETION))
SPLICE_SITE_ACCEPTOR=as.matrix(read.table(file=input.SPLICE_SITE_ACCEPTOR))
SPLICE_SITE_DONOR=as.matrix(read.table(file=input.SPLICE_SITE_DONOR))
UTR_5_PRIME=as.matrix(read.table(file=input.UTR_5_PRIME))
UTR_3_PRIME=as.matrix(read.table(file=input.UTR_3_PRIME))

out.EXON_DELETED=matrix(0,nrow=nrow(gene), ncol=2)
out.FRAME_SHIFT=matrix(0,nrow=nrow(gene), ncol=2)
out.CODON_CHANGE=matrix(0,nrow=nrow(gene), ncol=2)
out.UTR_5_DELETED=matrix(0,nrow=nrow(gene), ncol=2)
out.UTR_3_DELETED=matrix(0,nrow=nrow(gene), ncol=2)
out.CODON_INSERTION=matrix(0,nrow=nrow(gene), ncol=2)
out.CODON_CHANGE_PLUS_CODON_INSERTION=matrix(0,nrow=nrow(gene), ncol=2)
out.CODON_DELETION=matrix(0,nrow=nrow(gene), ncol=2)
out.CODON_CHANGE_PLUS_CODON_DELETION=matrix(0,nrow=nrow(gene), ncol=2)
out.SPLICE_SITE_ACCEPTOR=matrix(0,nrow=nrow(gene), ncol=2)
out.SPLICE_SITE_DONOR=matrix(0,nrow=nrow(gene), ncol=2)
out.UTR_5_PRIME=matrix(0,nrow=nrow(gene), ncol=2)
out.UTR_3_PRIME=matrix(0,nrow=nrow(gene), ncol=2)

out.EXON_DELETED[,1] = out.FRAME_SHIFT[,1] = out.CODON_CHANGE[,1] = out.UTR_5_DELETED[,1] = out.UTR_3_DELETED[,1] = out.CODON_INSERTION[,1] = out.CODON_CHANGE_PLUS_CODON_INSERTION[,1] = out.CODON_DELETION[,1] = out.CODON_CHANGE_PLUS_CODON_DELETION[,1] = out.SPLICE_SITE_ACCEPTOR[,1] = out.SPLICE_SITE_DONOR[,1] = out.UTR_5_PRIME[,1] = out.UTR_3_PRIME[,1] = as.character(gene[,1])

for (i in 1:nrow(EXON_DELETED))
{
	out.EXON_DELETED[out.EXON_DELETED[,1]==EXON_DELETED[i],2] = as.numeric(out.EXON_DELETED[out.EXON_DELETED[,1]==EXON_DELETED[i],2]) + 1
}
for (i in 1:nrow(FRAME_SHIFT))
{
	out.FRAME_SHIFT[out.FRAME_SHIFT[,1]==FRAME_SHIFT[i],2] = as.numeric(out.FRAME_SHIFT[out.FRAME_SHIFT[,1]==FRAME_SHIFT[i],2]) + 1
}
for (i in 1:nrow(CODON_CHANGE))
{
	out.CODON_CHANGE[out.CODON_CHANGE[,1]==CODON_CHANGE[i],2] = as.numeric(out.CODON_CHANGE[out.CODON_CHANGE[,1]==CODON_CHANGE[i],2]) + 1
}
for (i in 1:nrow(UTR_5_DELETED))
{
	out.UTR_5_DELETED[out.UTR_5_DELETED[,1]==UTR_5_DELETED[i],2] = as.numeric(out.UTR_5_DELETED[out.UTR_5_DELETED[,1]==UTR_5_DELETED[i],2]) + 1
}
for (i in 1:nrow(UTR_3_DELETED))
{
	out.UTR_3_DELETED[out.UTR_3_DELETED[,1]==UTR_3_DELETED[i],2] = as.numeric(out.UTR_3_DELETED[out.UTR_3_DELETED[,1]==UTR_3_DELETED[i],2]) + 1
}
for (i in 1:nrow(CODON_INSERTION))
{
	out.CODON_INSERTION[out.CODON_INSERTION[,1]==CODON_INSERTION[i],2] = as.numeric(out.CODON_INSERTION[out.CODON_INSERTION[,1]==CODON_INSERTION[i],2]) + 1
}
for (i in 1:nrow(CODON_CHANGE_PLUS_CODON_INSERTION))
{
	out.CODON_CHANGE_PLUS_CODON_INSERTION[out.CODON_CHANGE_PLUS_CODON_INSERTION[,1]==CODON_CHANGE_PLUS_CODON_INSERTION[i],2] = as.numeric(out.CODON_CHANGE_PLUS_CODON_INSERTION[out.CODON_CHANGE_PLUS_CODON_INSERTION[,1]==CODON_CHANGE_PLUS_CODON_INSERTION[i],2]) + 1
}
for (i in 1:nrow(CODON_DELETION))
{
	out.CODON_DELETION[out.CODON_DELETION[,1]==CODON_DELETION[i],2] = as.numeric(out.CODON_DELETION[out.CODON_DELETION[,1]==CODON_DELETION[i],2]) + 1
}
for (i in 1:nrow(CODON_CHANGE_PLUS_CODON_DELETION))
{
	out.CODON_CHANGE_PLUS_CODON_DELETION[out.CODON_CHANGE_PLUS_CODON_DELETION[,1]==CODON_CHANGE_PLUS_CODON_DELETION[i],2] = as.numeric(out.CODON_CHANGE_PLUS_CODON_DELETION[out.CODON_CHANGE_PLUS_CODON_DELETION[,1]==CODON_CHANGE_PLUS_CODON_DELETION[i],2]) + 1
}
for (i in 1:nrow(SPLICE_SITE_ACCEPTOR))
{
	out.SPLICE_SITE_ACCEPTOR[out.SPLICE_SITE_ACCEPTOR[,1]==SPLICE_SITE_ACCEPTOR[i],2] = as.numeric(out.SPLICE_SITE_ACCEPTOR[out.SPLICE_SITE_ACCEPTOR[,1]==SPLICE_SITE_ACCEPTOR[i],2]) + 1
}
for (i in 1:nrow(SPLICE_SITE_DONOR))
{
	out.SPLICE_SITE_DONOR[out.SPLICE_SITE_DONOR[,1]==SPLICE_SITE_DONOR[i],2] = as.numeric(out.SPLICE_SITE_DONOR[out.SPLICE_SITE_DONOR[,1]==SPLICE_SITE_DONOR[i],2]) + 1
}
for (i in 1:nrow(UTR_5_PRIME))
{
	out.UTR_5_PRIME[out.UTR_5_PRIME[,1]==UTR_5_PRIME[i],2] = as.numeric(out.UTR_5_PRIME[out.UTR_5_PRIME[,1]==UTR_5_PRIME[i],2]) + 1
}
for (i in 1:nrow(UTR_3_PRIME))
{
	out.UTR_3_PRIME[out.UTR_3_PRIME[,1]==UTR_3_PRIME[i],2] = as.numeric(out.UTR_3_PRIME[out.UTR_3_PRIME[,1]==UTR_3_PRIME[i],2]) + 1
}

write.table(out.EXON_DELETED, file = output.EXON_DELETED, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.FRAME_SHIFT, file = output.FRAME_SHIFT, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.CODON_CHANGE, file = output.CODON_CHANGE, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.UTR_5_DELETED, file = output.UTR_5_DELETED, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.UTR_3_DELETED, file = output.UTR_3_DELETED, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.CODON_INSERTION, file = output.CODON_INSERTION, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.CODON_CHANGE_PLUS_CODON_INSERTION, file = output.CODON_CHANGE_PLUS_CODON_INSERTION, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.CODON_DELETION, file = output.CODON_DELETION, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.CODON_CHANGE_PLUS_CODON_DELETION, file = output.CODON_CHANGE_PLUS_CODON_DELETION, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.SPLICE_SITE_ACCEPTOR, file = output.SPLICE_SITE_ACCEPTOR, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.SPLICE_SITE_DONOR, file = output.SPLICE_SITE_DONOR, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.UTR_5_PRIME, file = output.UTR_5_PRIME, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.UTR_3_PRIME, file = output.UTR_3_PRIME, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
