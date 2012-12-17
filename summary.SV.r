stdin = commandArgs() 

for (i in 1:1)
{
	if (length(stdin)==16)
	{
	input.gene=as.character(stdin[6])
	input.ITX=as.character(stdin[7])
	input.INV=as.character(stdin[8])
	input.DEL=as.character(stdin[9])
	input.INS=as.character(stdin[10])
	input.CTX=as.character(stdin[11])
	output.ITX=as.character(stdin[12])
	output.INV=as.character(stdin[13])
	output.DEL=as.character(stdin[14])
	output.INS=as.character(stdin[15])
	output.CTX=as.character(stdin[16])
	}
}
gene=as.matrix(read.table(file=input.gene))
ITX=as.matrix(read.table(file=input.ITX))
INV=as.matrix(read.table(file=input.INV))
DEL=as.matrix(read.table(file=input.DEL))
INS=as.matrix(read.table(file=input.INS))
CTX=as.matrix(read.table(file=input.CTX))
out.ITX=matrix(0,nrow=nrow(gene), ncol=2)
out.INV=matrix(0,nrow=nrow(gene), ncol=2)
out.DEL=matrix(0,nrow=nrow(gene), ncol=2)
out.INS=matrix(0,nrow=nrow(gene), ncol=2)
out.CTX=matrix(0,nrow=nrow(gene), ncol=2)

out.ITX[,1] = out.INV[,1] = out.DEL[,1] = out.INS[,1] = out.CTX[,1] = as.character(gene[,1])

for (i in 1:nrow(ITX))
{
	out.ITX[out.ITX[,1]==ITX[i,1],2] = as.numeric(out.ITX[out.ITX[,1]==ITX[i,1],2]) + 1
}
for (i in 1:nrow(INV))
{
	out.INV[out.INV[,1]==INV[i,1],2] = as.numeric(out.INV[out.INV[,1]==INV[i,1],2]) + 1
}
for (i in 1:nrow(DEL))
{
	out.DEL[out.DEL[,1]==DEL[i,1],2] = as.numeric(out.DEL[out.DEL[,1]==DEL[i,1],2]) + 1
}
for (i in 1:nrow(INS))
{
	out.INS[out.INS[,1]==INS[i,1],2] = as.numeric(out.INS[out.INS[,1]==INS[i,1],2]) + 1
}
for (i in 1:nrow(CTX))
{
	out.CTX[out.CTX[,1]==CTX[i,1],2] = as.numeric(out.CTX[out.CTX[,1]==CTX[i,1],2]) + 1
}

# for (i in 1:nrow(ITX))
# {
	# if (as.character(ITX[i,2]) == as.character(ITX[i,3]))
	# {
		# out.ITX[out.ITX[,1]==ITX[i,2],2] = as.numeric(out.ITX[out.ITX[,1]==ITX[i,2],2]) + 1
	# }
	# else if (as.character(ITX[i,2]) == "NOGENE")
	# {
		# out.ITX[out.ITX[,1]==ITX[i,3],2] = as.numeric(out.ITX[out.ITX[,1]==ITX[i,3],2]) + 1
	# }
	# else if (as.character(ITX[i,3]) == "NOGENE")
	# {
		# out.ITX[out.ITX[,1]==ITX[i,2],2] = as.numeric(out.ITX[out.ITX[,1]==ITX[i,2],2]) + 1
	# }
	# else
	# {
		# out.ITX[out.ITX[,1]==ITX[i,2],2] = as.numeric(out.ITX[out.ITX[,1]==ITX[i,2],2]) + 1
		# out.ITX[out.ITX[,1]==ITX[i,3],2] = as.numeric(out.ITX[out.ITX[,1]==ITX[i,3],2]) + 1
	# }
# }

# for (i in 1:nrow(INV))
# {
	# if (as.character(INV[i,2]) == as.character(INV[i,3]))
	# {
		# out.INV[out.INV[,1]==INV[i,2],2] = as.numeric(out.INV[out.INV[,1]==INV[i,2],2]) + 1
	# }
	# else if (as.character(INV[i,2]) == "NOGENE")
	# {
		# out.INV[out.INV[,1]==INV[i,3],2] = as.numeric(out.INV[out.INV[,1]==INV[i,3],2]) + 1
	# }
	# else if (as.character(INV[i,3]) == "NOGENE")
	# {
		# out.INV[out.INV[,1]==INV[i,2],2] = as.numeric(out.INV[out.INV[,1]==INV[i,2],2]) + 1
	# }
	# else
	# {
		# out.INV[out.INV[,1]==INV[i,2],2] = as.numeric(out.INV[out.INV[,1]==INV[i,2],2]) + 1
		# out.INV[out.INV[,1]==INV[i,3],2] = as.numeric(out.INV[out.INV[,1]==INV[i,3],2]) + 1
	# }
# }

# for (i in 1:nrow(DEL))
# {
	# if (as.character(DEL[i,2]) == as.character(DEL[i,3]))
	# {
		# out.DEL[out.DEL[,1]==DEL[i,2],2] = as.numeric(out.DEL[out.DEL[,1]==DEL[i,2],2]) + 1
	# }
	# else if (as.character(DEL[i,2]) == "NOGENE")
	# {
		# out.DEL[out.DEL[,1]==DEL[i,3],2] = as.numeric(out.DEL[out.DEL[,1]==DEL[i,3],2]) + 1
	# }
	# else if (as.character(DEL[i,3]) == "NOGENE")
	# {
		# out.DEL[out.DEL[,1]==DEL[i,2],2] = as.numeric(out.DEL[out.DEL[,1]==DEL[i,2],2]) + 1
	# }
	# else
	# {
		# out.DEL[out.DEL[,1]==DEL[i,2],2] = as.numeric(out.DEL[out.DEL[,1]==DEL[i,2],2]) + 1
		# out.DEL[out.DEL[,1]==DEL[i,3],2] = as.numeric(out.DEL[out.DEL[,1]==DEL[i,3],2]) + 1
	# }
# }

# for (i in 1:nrow(INS))
# {
	# if (as.character(INS[i,2]) == as.character(INS[i,3]))
	# {
		# out.INS[out.INS[,1]==INS[i,2],2] = as.numeric(out.INS[out.INS[,1]==INS[i,2],2]) + 1
	# }
	# else if (as.character(INS[i,2]) == "NOGENE")
	# {
		# out.INS[out.INS[,1]==INS[i,3],2] = as.numeric(out.INS[out.INS[,1]==INS[i,3],2]) + 1
	# }
	# else if (as.character(INS[i,3]) == "NOGENE")
	# {
		# out.INS[out.INS[,1]==INS[i,2],2] = as.numeric(out.INS[out.INS[,1]==INS[i,2],2]) + 1
	# }
	# else
	# {
		# out.INS[out.INS[,1]==INS[i,2],2] = as.numeric(out.INS[out.INS[,1]==INS[i,2],2]) + 1
		# out.INS[out.INS[,1]==INS[i,3],2] = as.numeric(out.INS[out.INS[,1]==INS[i,3],2]) + 1
	# }
# }

# for (i in 1:nrow(CTX))
# {
	# if (as.character(CTX[i,2]) == as.character(CTX[i,3]))
	# {
		# out.CTX[out.CTX[,1]==CTX[i,2],2] = as.numeric(out.CTX[out.CTX[,1]==CTX[i,2],2]) + 1
	# }
	# else if (as.character(CTX[i,2]) == "NOGENE")
	# {
		# out.CTX[out.CTX[,1]==CTX[i,3],2] = as.numeric(out.CTX[out.CTX[,1]==CTX[i,3],2]) + 1
	# }
	# else if (as.character(CTX[i,3]) == "NOGENE")
	# {
		# out.CTX[out.CTX[,1]==CTX[i,2],2] = as.numeric(out.CTX[out.CTX[,1]==CTX[i,2],2]) + 1
	# }
	# else
	# {
		# out.CTX[out.CTX[,1]==CTX[i,2],2] = as.numeric(out.CTX[out.CTX[,1]==CTX[i,2],2]) + 1
		# out.CTX[out.CTX[,1]==CTX[i,3],2] = as.numeric(out.CTX[out.CTX[,1]==CTX[i,3],2]) + 1
	# }
# }

write.table(out.ITX, file = output.ITX, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.INV, file = output.INV, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.DEL, file = output.DEL, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.INS, file = output.INS, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(out.CTX, file = output.CTX, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
