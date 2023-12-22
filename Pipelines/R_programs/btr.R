#Load up the library
library(BTR)

#Sets the directory
setwd("/Users/mbrown/Desktop/Research/BTR")

set.seed(0)

#Reads in the matrix
mat <- read.table('Test2.txt', header = FALSE, sep = '\t')
mat <- initialise_raw_data(mat, max_expr = "low")
mat <- t(mat)
#Gives it headings
colnames(mat) <- paste("Gene", seq(1:ncol(mat)), sep = "_")
rownames(mat) <- paste("Cell", seq(1:nrow(mat)), sep = "_")

#Filters out genes and cells
gene_ind = c("Gene_1", "Gene_2", "Gene_3", "Gene_4", "Gene_5", "Gene_6")
mat = mat[, gene_ind]
cell_ind = grepl("cmp", rownames(mat)) | grepl("gmp", rownames(mat)) | grepl("mep", 
                                                                                 rownames(mat))

mat = as.data.frame(mat)

FM = model_train(mat, max_varperrule = 2, verbose = T)
plotBM(FM)
