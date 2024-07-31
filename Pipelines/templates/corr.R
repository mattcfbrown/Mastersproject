#I am now going to try and use Pearson correlations

#The library
library(Hmisc)
library(MASS)

args = commandArgs(trailingOnly = TRUE)
infile = args[1]
id = args[2]

#The data:
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

output_corr <- rcorr(trans_mat, type="pearson")
output_corr <- output_corr$r

#Diagonal convert to 0 (as no internal regulation)
for (k in 1:length(gene_names)){
  output_corr[k,k] = 0
}

output_corr <- abs(output_corr)
sorted_output <- sort(unlist(as.list(output_corr)),decreasing = TRUE)
threshold <- sorted_output[round(0.05*length(sorted_output))]

#Save the weights:
rownames(output_corr) <- NULL
colnames(output_corr) <- NULL
name <- paste0("corr_weight_",as.character(id))
write.matrix(output_corr,file=paste0(name,".csv"), sep = ',')

guess <- matrix(0,length(gene_names),length(gene_names))
#Now get the guess
for (a in 1:length(gene_names)){
  for (b in 1:length(gene_names)){
    if (output_corr[a,b] >= threshold){
      guess[a,b] = 1
    }
  }
}

#Save the guess
rownames(guess) <- NULL
colnames(guess) <- NULL
name <- paste0("corr_guess_",as.character(id))
write.matrix(guess,file=paste0(name,".csv"), sep = ',')
