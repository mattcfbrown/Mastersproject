#Here I will be testing the capabilities of the TIGRESS R package

#Library needed
library(tigress)
library(MASS)

args = commandArgs(trailingOnly = TRUE)
infile = args[1]
id = args[2]

#This here is an example, once I have it working I will operate on other data
mat <- as.matrix(read.table(infile,header=FALSE,sep="\t"))
trans_mat <- t(mat)

#Get the gene and cell names required
gene_names <- c()
for (i in 1:ncol(trans_mat)){
  name <- paste0("T",as.character(i))
  gene_names <- c(gene_names,name)
}
cell_names <- c()
for (j in 1:nrow(trans_mat)){
  name <- paste0("Cell",as.character(j))
  cell_names <- c(cell_names,name)
}

colnames(trans_mat) <- gene_names
rownames(trans_mat) <- cell_names

#Now try to run Tigress
output <- tigress(trans_mat, allsteps = FALSE)

#We firstly save the weights
rownames(output) <- NULL
colnames(output) <- NULL
name <- paste0("tigress_weight_",as.character(id))
write.matrix(output,file=paste0(name,".csv"), sep = ',')

sorted_output <- sort(unlist(as.list(output)),decreasing = TRUE)
threshold <- sorted_output[round(0.05*length(sorted_output))]

#We now save the guess
guess <- matrix(0,length(gene_names),length(gene_names))
#Now get the guess
for (a in 1:length(gene_names)){
  for (b in 1:length(gene_names)){
    if (output[a,b] >= threshold){
      guess[a,b] = 1
    }
  }
}

#Save the guess
rownames(guess) <- NULL
colnames(guess) <- NULL
name <- paste0("tigress_guess_",as.character(id))
write.matrix(guess,file=paste0(name,".csv"), sep = ',')


