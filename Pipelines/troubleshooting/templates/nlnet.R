#Load up the library
library(nlnet)

#Sets the directory
#setwd("/Users/mbrown/Desktop/Research/NLNET")

#Now get an example matrix
x <- data.gen(n.genes = 30, n.samples = 400)

nlnet(x[["data"]])

# pass this as a input file to the script and inject below rather than hard code
mat <- read.table('Test2.txt', header = FALSE, sep = '\t')

data_matrix = data.matrix(mat)

nlnet(data_matrix)
