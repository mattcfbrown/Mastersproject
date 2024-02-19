#Code used from
#https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html#using-slingshot

#Library to use
library(slingshot)
library(mclust, quietly = TRUE)

#Read in the count matrix
args = commandArgs(trailingOnly = TRUE)

counts <- read.csv(args[1], header=FALSE, sep= ',') 
counts = as.matrix(counts)
num_genes <- dim(counts)[1]
num_cells <- dim(counts)[2]

#Converts it into a form we can use
rownames(counts) <- paste0('Genes',1:num_genes)
colnames(counts) <- paste0('Cell',1:num_cells)
toy <- SingleCellExperiment(assays = List(counts = counts))


#Normalisation
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(toy)$norm <- FQnorm(assays(toy)$counts)

#Perform PCA on it
pca <- prcomp(t(log1p(assays(toy)$norm)), scale. = FALSE)
reduced <- pca$x[,1:2]
reducedDims(toy) <- SimpleList(PCA = reduced)


#Cell clustering time
GMM_cluster <- Mclust(reduced)$classification
colData(toy)$GMM <- GMM_cluster

#Gets pseudotime
toy <- slingshot(toy, clusterLabels = 'GMM', reducedDim = 'PCA')
pseudotime <- matrix(slingAvgPseudotime(toy))
pseudotime <- cbind(pseudotime, GMM_cluster)

#Now normalise the data between 0 to 1
min <- min(pseudotime[,1])
max <- max(pseudotime[,1])
x <- dim(pseudotime)[1]
norm <- rep(0,x)

#https://www.statology.org/normalize-data-between-0-and-1/
for (i in 1:x) {
  norm[i] <- (pseudotime[i,1]-min)/(max-min)
}

pseudotime <- cbind(pseudotime,norm)
pseudotime <- pseudotime[order(pseudotime[,2]),]

write.table(pseudotime[,2:3], "mytable_R.txt", col.name = FALSE, row.name = FALSE, sep = "\t")
