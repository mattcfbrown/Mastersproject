#Load up the library
library(GENIE3)
set.seed(123)

args = commandArgs(trailingOnly = TRUE)

#Reads in the data
mat <- read.table(args[1], header = FALSE, sep = '\t')
num_gene <- nrow(mat)

#Gives the names needed
rownames(mat) <- paste("Gene", 1:num_gene, sep="")
colnames(mat) <- paste("Sample", 1:ncol(mat), sep="")

#Converts mat into a data matrix
mat <- data.matrix(mat)

#Performs the GENIE3 algorithm
weightMat <- GENIE3(mat)
linklist <- getLinkList(weightMat)
linkList_threshold <- getLinkList(weightMat, threshold=args[2])

print(linklist_threshold)