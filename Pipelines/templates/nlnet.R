#Load up the library
library(nlnet)

args = commandArgs(trailingOnly = TRUE)


#Now get an example matrix
#x <- data.gen(n.genes = 30, n.samples = 400)

#nlnet(x[["data"]])

# pass this as a input file to the script and inject below rather than hard code
mat <- read.table(args, header = FALSE, sep = '\t')

data_matrix = data.matrix(mat)

test <- nlnet(data_matrix)

num_genes = length(test$graph)
values <- c()
for (x in 1:num_genes){
  values <- append(values,test$graph[x])
}

adjacency <- matrix(values,nrow=num_genes,ncol = num_genes)
adjacency