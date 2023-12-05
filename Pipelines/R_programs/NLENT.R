#Load up the library
library(nlnet)

#Sets the directory
setwd("/Users/mbrown/Desktop/Research/NLNET")

#Now get an example matrix
x <- data.gen(n.genes = 30, n.samples = 400)

nlnet(x[["data"]])

mat <- read.table('Test2.txt', header = FALSE, sep = '\t')

mat_convert = data.matrix(mat)

nlnet(mat_convert)
