stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==7)
	{
	input=as.character(stdin[6])
	output=as.character(stdin[7])
	}
}
list=as.matrix(read.table(file=input))
for (i in nrow(list):1)
{
	if (list[i,5]=="NOGENE")
	{
		if (list[i,10]=="NOGENE")
		{
			list <- list[-i,]
		}
	}
}
write.table(list, file = output, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
