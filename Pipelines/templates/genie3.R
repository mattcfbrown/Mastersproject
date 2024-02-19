#Load up the library
library(GENIE3)
set.seed(123)

# assign variables from command line
args = commandArgs(trailingOnly = TRUE)
infile = args[1]
threshold = as.numeric(args[2])

#Reads in the data
mat <- read.table(infile, header = FALSE, sep = '\t')
num_gene <- nrow(mat)

#Gives the names needed
rownames(mat) <- paste("Gene", 1:num_gene, sep="")
colnames(mat) <- paste("Sample", 1:ncol(mat), sep="")

#Converts mat into a data matrix
mat <- data.matrix(mat)

#Performs the GENIE3 algorithm
weightMat <- GENIE3(mat)
linklist <- getLinkList(weightMat)
linkList_threshold <- getLinkList(weightMat, threshold=threshold)

print(linkList_threshold)