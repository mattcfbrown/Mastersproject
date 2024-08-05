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
output_spear <- rcorr(trans_mat, type="spearman")
output_corr <- output_corr$r
output_spear <-output_spear$r

#Diagonal convert to 0 (as no internal regulation)
for (k in 1:length(gene_names)){
  output_corr[k,k] = 0
  output_spear[k,k] = 0  
}

output_corr <- abs(output_corr)
output_spear <- abs(output_spear)
sorted_output <- sort(unlist(as.list(output_corr)),decreasing = TRUE)
sorted_spear <- sort(unlist(as.list(output_spear)),decreasing = TRUE)
threshold <- sorted_output[round(0.05*length(sorted_output))]
thresh_spear <- sorted_spear[round(0.05*length(sorted_spear))]

#Min-max normalise the data
max <- max(output_corr)
min <- min(output_corr)
norm_matrix <- matrix(0,length(gene_names),length(gene_names))
for (i in 1:length(gene_names)){
  for (j in 1:length(gene_names)){
    norm_matrix[i,j] = (output_corr[i,j]-min)/(max-min)
  }
}

max <- max(output_spear)
min <- min(output_spear)
spear_matrix <- matrix(0,length(gene_names),length(gene_names))
for (i in 1:length(gene_names)){
  for (j in 1:length(gene_names)){
    spear_matrix[i,j] = (output_spear[i,j]-min)/(max-min)
  }
}

#Save the weights:
rownames(norm_matrix) <- NULL
colnames(norm_matrix) <- NULL
rownames(spear_matrix) <- NULL
colnames(spear_matrix) <- NULL
name <- paste0("pearson_weight_",as.character(id))
write.matrix(norm_matrix,file=paste0(name,".csv"), sep = ',')
name <- paste0("spear_weight_",as.character(id))
write.matrix(spear_matrix,file=paste0(name,".csv"), sep = ',')

#--------------------Now the guesses-----------------------
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
name <- paste0("pearson_guess_",as.character(id))
write.matrix(guess,file=paste0(name,".csv"), sep = ',')

#Spearman
guess_spear <- matrix(0,length(gene_names),length(gene_names))
#Now get the guess
for (a in 1:length(gene_names)){
  for (b in 1:length(gene_names)){
    if (output_spear[a,b] >= thresh_spear){
      guess_spear[a,b] = 1
    }
  }
}

#Save the guess
rownames(guess_spear) <- NULL
colnames(guess_spear) <- NULL
name <- paste0("spear_guess_",as.character(id))
write.matrix(guess_spear,file=paste0(name,".csv"), sep = ',')
