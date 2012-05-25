stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==36)
	{
	input.gene=as.character(stdin[6])
	input.SPLICE_SITE_ACCEPTOR=as.character(stdin[7])
	input.SPLICE_SITE_DONOR=as.character(stdin[8])
	input.START_LOST=as.character(stdin[9])
	input.STOP_GAINED=as.character(stdin[10])
	input.STOP_LOST=as.character(stdin[11])
	input.RARE_AMINO_ACID=as.character(stdin[12])
	input.NON_SYNONYMOUS_CODING=as.character(stdin[13])
	input.SYNONYMOUS_START=as.character(stdin[14])
	input.NON_SYNONYMOUS_START=as.character(stdin[15])
	input.START_GAINED=as.character(stdin[16])
	input.SYNONYMOUS_CODING=as.character(stdin[17])
	input.SYNONYMOUS_STOP=as.character(stdin[18])
	input.NON_SYNONYMOUS_STOP=as.character(stdin[19])
	input.UTR_5_PRIME=as.character(stdin[20])
	input.UTR_3_PRIME=as.character(stdin[21])
	output.SPLICE_SITE_ACCEPTOR=as.character(stdin[22])
	output.SPLICE_SITE_DONOR=as.character(stdin[23])
	output.START_LOST=as.character(stdin[24])
	output.STOP_GAINED=as.character(stdin[25])
	output.STOP_LOST=as.character(stdin[26])
	output.RARE_AMINO_ACID=as.character(stdin[27])
	output.NON_SYNONYMOUS_CODING=as.character(stdin[28])
	output.SYNONYMOUS_START=as.character(stdin[29])
	output.NON_SYNONYMOUS_START=as.character(stdin[30])
	output.START_GAINED=as.character(stdin[31])
	output.SYNONYMOUS_CODING=as.character(stdin[32])
	output.SYNONYMOUS_STOP=as.character(stdin[33])
	output.NON_SYNONYMOUS_STOP=as.character(stdin[34])
	output.UTR_5_PRIME=as.character(stdin[35])
	output.UTR_3_PRIME=as.character(stdin[36])
	}
}
gene=as.matrix(read.table(file=input.gene))
SPLICE_SITE_ACCEPTOR=as.matrix(read.table(file=input.SPLICE_SITE_ACCEPTOR))
SPLICE_SITE_DONOR=as.matrix(read.table(file=input.SPLICE_SITE_DONOR))
START_LOST=as.matrix(read.table(file=input.START_LOST))
STOP_GAINED=as.matrix(read.table(file=input.STOP_GAINED))
STOP_LOST=as.matrix(read.table(file=input.STOP_LOST))
RARE_AMINO_ACID=as.matrix(read.table(file=input.RARE_AMINO_ACID))
NON_SYNONYMOUS_CODING=as.matrix(read.table(file=input.NON_SYNONYMOUS_CODING))
SYNONYMOUS_START=as.matrix(read.table(file=input.SYNONYMOUS_START))
NON_SYNONYMOUS_START=as.matrix(read.table(file=input.NON_SYNONYMOUS_START))
START_GAINED=as.matrix(read.table(file=input.START_GAINED))
SYNONYMOUS_CODING=as.matrix(read.table(file=input.SYNONYMOUS_CODING))
SYNONYMOUS_STOP=as.matrix(read.table(file=input.SYNONYMOUS_STOP))
NON_SYNONYMOUS_STOP=as.matrix(read.table(file=input.NON_SYNONYMOUS_STOP))
UTR_5_PRIME=as.matrix(read.table(file=input.UTR_5_PRIME))
UTR_3_PRIME=as.matrix(read.table(file=input.UTR_3_PRIME))
out.SPLICE_SITE_ACCEPTOR=matrix(0,nrow=nrow(gene), ncol=2)
out.SPLICE_SITE_DONOR=matrix(0,nrow=nrow(gene), ncol=2)
out.START_LOST=matrix(0,nrow=nrow(gene), ncol=2)
out.STOP_GAINED=matrix(0,nrow=nrow(gene), ncol=2)
out.STOP_LOST=matrix(0,nrow=nrow(gene), ncol=2)
out.RARE_AMINO_ACID=matrix(0,nrow=nrow(gene), ncol=2)
out.NON_SYNONYMOUS_CODING=matrix(0,nrow=nrow(gene), ncol=2)
out.SYNONYMOUS_START=matrix(0,nrow=nrow(gene), ncol=2)
out.NON_SYNONYMOUS_START=matrix(0,nrow=nrow(gene), ncol=2)
out.START_GAINED=matrix(0,nrow=nrow(gene), ncol=2)
out.SYNONYMOUS_CODING=matrix(0,nrow=nrow(gene), ncol=2)
out.SYNONYMOUS_STOP=matrix(0,nrow=nrow(gene), ncol=2)
out.NON_SYNONYMOUS_STOP=matrix(0,nrow=nrow(gene), ncol=2)
out.UTR_5_PRIME=matrix(0,nrow=nrow(gene), ncol=2)
out.UTR_3_PRIME=matrix(0,nrow=nrow(gene), ncol=2)

out.SPLICE_SITE_ACCEPTOR[,1] = out.SPLICE_SITE_DONOR[,1] = out.START_LOST[,1] = out.STOP_GAINED[,1] = out.STOP_LOST[,1] = out.RARE_AMINO_ACID[,1] = out.NON_SYNONYMOUS_CODING[,1] = out.SYNONYMOUS_START[,1] =out.NON_SYNONYMOUS_START[,1] = out.START_GAINED[,1] = out.SYNONYMOUS_CODING[,1] = out.SYNONYMOUS_STOP[,1] = out.NON_SYNONYMOUS_STOP[,1] = out.UTR_5_PRIME[,1] = out.UTR_3_PRIME[,1] = as.character(gene[,1])

for (i in 1:nrow(SPLICE_SITE_ACCEPTOR))
{
	out.SPLICE_SITE_ACCEPTOR[out.SPLICE_SITE_ACCEPTOR[,1]==SPLICE_SITE_ACCEPTOR[i],2] = as.numeric(out.SPLICE_SITE_ACCEPTOR[out.SPLICE_SITE_ACCEPTOR[,1]==SPLICE_SITE_ACCEPTOR[i],2]) + 1
}
for (i in 1:nrow(SPLICE_SITE_DONOR))
{
	out.SPLICE_SITE_DONOR[out.SPLICE_SITE_DONOR[,1]==SPLICE_SITE_DONOR[i],2] = as.numeric(out.SPLICE_SITE_DONOR[out.SPLICE_SITE_DONOR[,1]==SPLICE_SITE_DONOR[i],2]) + 1
}
for (i in 1:nrow(START_LOST))
{
	out.START_LOST[out.START_LOST[,1]==START_LOST[i],2] = as.numeric(out.START_LOST[out.START_LOST[,1]==START_LOST[i],2]) + 1
}
for (i in 1:nrow(STOP_GAINED))
{
	out.STOP_GAINED[out.STOP_GAINED[,1]==STOP_GAINED[i],2] = as.numeric(out.STOP_GAINED[out.STOP_GAINED[,1]==STOP_GAINED[i],2]) + 1
}

for (i in 1:nrow(STOP_LOST))
{
	out.STOP_LOST[out.STOP_LOST[,1]==STOP_LOST[i],2] = as.numeric(out.STOP_LOST[out.STOP_LOST[,1]==STOP_LOST[i],2]) + 1
}
for (i in 1:nrow(RARE_AMINO_ACID))
{
	out.RARE_AMINO_ACID[out.RARE_AMINO_ACID[,1]==RARE_AMINO_ACID[i],2] = as.numeric(out.RARE_AMINO_ACID[out.RARE_AMINO_ACID[,1]==RARE_AMINO_ACID[i],2]) + 1
}
for (i in 1:nrow(NON_SYNONYMOUS_CODING))
{
	out.NON_SYNONYMOUS_CODING[out.NON_SYNONYMOUS_CODING[,1]==NON_SYNONYMOUS_CODING[i],2] = as.numeric(out.NON_SYNONYMOUS_CODING[out.NON_SYNONYMOUS_CODING[,1]==NON_SYNONYMOUS_CODING[i],2]) + 1
}
for (i in 1:nrow(SYNONYMOUS_START))
{
	out.SYNONYMOUS_START[out.SYNONYMOUS_START[,1]==SYNONYMOUS_START[i],2] = as.numeric(out.SYNONYMOUS_START[out.SYNONYMOUS_START[,1]==SYNONYMOUS_START[i],2]) + 1
}
for (i in 1:nrow(NON_SYNONYMOUS_START))
{
	out.NON_SYNONYMOUS_START[out.NON_SYNONYMOUS_START[,1]==NON_SYNONYMOUS_START[i],2] = as.numeric(out.NON_SYNONYMOUS_START[out.NON_SYNONYMOUS_START[,1]==NON_SYNONYMOUS_START[i],2]) + 1
}
for (i in 1:nrow(START_GAINED))
{
	out.START_GAINED[out.START_GAINED[,1]==START_GAINED[i],2] = as.numeric(out.START_GAINED[out.START_GAINED[,1]==START_GAINED[i],2]) + 1
}
for (i in 1:nrow(SYNONYMOUS_CODING))
{
	out.SYNONYMOUS_CODING[out.SYNONYMOUS_CODING[,1]==SYNONYMOUS_CODING[i],2] = as.numeric(out.SYNONYMOUS_CODING[out.SYNONYMOUS_CODING[,1]==SYNONYMOUS_CODING[i],2]) + 1
}
for (i in 1:nrow(SYNONYMOUS_STOP))
{
	out.SYNONYMOUS_STOP[out.SYNONYMOUS_STOP[,1]==SYNONYMOUS_STOP[i],2] = as.numeric(out.SYNONYMOUS_STOP[out.SYNONYMOUS_STOP[,1]==SYNONYMOUS_STOP[i],2]) + 1
}
for (i in 1:nrow(NON_SYNONYMOUS_STOP))
{
	out.NON_SYNONYMOUS_STOP[out.NON_SYNONYMOUS_STOP[,1]==NON_SYNONYMOUS_STOP[i],2] = as.numeric(out.NON_SYNONYMOUS_STOP[out.NON_SYNONYMOUS_STOP[,1]==NON_SYNONYMOUS_STOP[i],2]) + 1
}
for (i in 1:nrow(UTR_5_PRIME))
{
	out.UTR_5_PRIME[out.UTR_5_PRIME[,1]==UTR_5_PRIME[i],2] = as.numeric(out.UTR_5_PRIME[out.UTR_5_PRIME[,1]==UTR_5_PRIME[i],2]) + 1
}
for (i in 1:nrow(UTR_3_PRIME))
{
	out.UTR_3_PRIME[out.UTR_3_PRIME[,1]==UTR_3_PRIME[i],2] = as.numeric(out.UTR_3_PRIME[out.UTR_3_PRIME[,1]==UTR_3_PRIME[i],2]) + 1
}

write.table(out.SPLICE_SITE_ACCEPTOR, file = output.SPLICE_SITE_ACCEPTOR, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.SPLICE_SITE_DONOR, file = output.SPLICE_SITE_DONOR, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.START_LOST, file = output.START_LOST, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.STOP_GAINED, file = output.STOP_GAINED, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.STOP_LOST, file = output.STOP_LOST, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.RARE_AMINO_ACID, file = output.RARE_AMINO_ACID, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.NON_SYNONYMOUS_CODING, file = output.NON_SYNONYMOUS_CODING, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.SYNONYMOUS_START, file = output.SYNONYMOUS_START, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.SYNONYMOUS_START, file = output.SYNONYMOUS_START, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.NON_SYNONYMOUS_START, file = output.NON_SYNONYMOUS_START, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.START_GAINED, file = output.START_GAINED, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.SYNONYMOUS_CODING, file = output.SYNONYMOUS_CODING, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.SYNONYMOUS_STOP, file = output.SYNONYMOUS_STOP, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.UTR_5_PRIME, file = output.UTR_5_PRIME, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.NON_SYNONYMOUS_STOP, file = output.NON_SYNONYMOUS_STOP, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.UTR_3_PRIME, file = output.UTR_3_PRIME, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
