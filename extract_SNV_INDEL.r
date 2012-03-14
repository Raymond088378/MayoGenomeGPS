stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==8)
	{
	input=as.character(stdin[6])
	output.INDEL=as.character(stdin[7])
	output.SNV=as.character(stdin[8])
	}
}
list=read.table(file=input)
numrows=nrow(list)
INDEL=matrix(NA, nrow=nrow(list), ncol=ncol(list))
SNV=matrix(NA, nrow=nrow(list), ncol=ncol(list))

for (i in 1:numrows)
{
if (nchar(as.character(list[i,5])) > 1)
{
INDEL[i,1]=as.character(list[i,1])
INDEL[i,2]=as.character(list[i,2])
INDEL[i,3]=as.character(list[i,3])
INDEL[i,4]=as.character(list[i,4])
INDEL[i,5]=as.character(list[i,5])
INDEL[i,6]=as.character(list[i,6])
INDEL[i,7]=as.character(list[i,7])
INDEL[i,8]=as.character(list[i,8])
INDEL[i,9]=as.character(list[i,9])
INDEL[i,10]=as.character(list[i,10])
INDEL[i,11]=as.character(list[i,11])
}
else if (nchar(as.character(list[i,6])) > 1)
{
INDEL[i,1]=as.character(list[i,1])
INDEL[i,2]=as.character(list[i,2])
INDEL[i,3]=as.character(list[i,3])
INDEL[i,4]=as.character(list[i,4])
INDEL[i,5]=as.character(list[i,5])
INDEL[i,6]=as.character(list[i,6])
INDEL[i,7]=as.character(list[i,7])
INDEL[i,8]=as.character(list[i,8])
INDEL[i,9]=as.character(list[i,9])
INDEL[i,10]=as.character(list[i,10])
INDEL[i,11]=as.character(list[i,11])
}
}

for (i in 1:numrows)
{
if ((nchar(as.character(list[i,5])) == 1) && (nchar(as.character(list[i,6])) == 1))
{
SNV[i,1]=as.character(list[i,1])
SNV[i,2]=as.character(list[i,2])
SNV[i,3]=as.character(list[i,3])
SNV[i,4]=as.character(list[i,4])
SNV[i,5]=as.character(list[i,5])
SNV[i,6]=as.character(list[i,6])
SNV[i,7]=as.character(list[i,7])
SNV[i,8]=as.character(list[i,8])
SNV[i,9]=as.character(list[i,9])
SNV[i,10]=as.character(list[i,10])
SNV[i,11]=as.character(list[i,11])
}
}
INDEL <- INDEL[complete.cases(INDEL),]
SNV <- SNV[complete.cases(SNV),]

write.table(INDEL, file = output.INDEL, append = FALSE, quote = FALSE, sep = "*", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(SNV, file = output.SNV, append = FALSE, quote = FALSE, sep = "*", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)