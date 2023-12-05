#Load up the library
library(GENIE3)

set.seed(123)

#Sets the directory
setwd("/Users/mbrown/Desktop/Research/GENIE3")

mat <- read.matrix('Test3.txt', header = FALSE, sep = '\t')
mat_df<-data.frame(mat)
apply(as.matrix.noquote(mat_df),2,as.numeric)

rownames(mat) <- paste("Gene", 1:25, sep="")
colnames(mat) <- paste("Sample", 1:2000, sep="")
head(mat_df)

mat <- data.matrix(mat)
class(mat)

weightMat <- GENIE3(mat)

linklist <- getLinkList(weightMat)

linkList_threshold <- getLinkList(weightMat, threshold=0.1)
