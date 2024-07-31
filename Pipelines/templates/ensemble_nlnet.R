#This script takes part of the nlnet function and tries to output the FDR matrix

normrow <- function(array) {
  m <- apply(array, 1, mean, na.rm = T)
  s <- apply(array, 1, sd, na.rm = T)
  array <- (array - m)/s
  return(array)
}

gene.specific.null <- function(array, B = 500) {
  null.mat <- matrix(0, nrow = nrow(array), ncol = B)
  l <- ncol(array)
  d.array <- array[, 1:(l - 1)]
  for (i in 1:B) {
    this.order <- sample(l, l, replace = FALSE)
    for (j in 1:(l - 1)) d.array[, j] <- abs(array[, 
                                                   this.order[j + 1]] - array[, this.order[j]])
    null.mat[, i] <- apply(d.array, 1, sum)
  }
  r <- cbind(apply(null.mat, 1, mean), apply(null.mat, 
                                             1, sd))
  return(r)
}

scol.matrix.order <- function(array, x) {
  if (is.null(nrow(array)) | nrow(array) == 1) {
    array <- as.vector(array)
    array <- array[order(x)]
    d <- array[2:length(array)] - array[1:(length(array) - 
                                             1)]
    dd <- sum(abs(d), na.rm = T)
  }
  else {
    array <- array[, order(x)]
    d <- array[, 2:ncol(array)] - array[, 1:(ncol(array) - 
                                               1)]
    dd <- apply(abs(d), 1, sum, na.rm = T)
  }
  return(dd)
}

scol.matrix <- function(a, direction = 2) {
  rdmat <- matrix(0, ncol = nrow(a), nrow = nrow(a))
  for (j in 1:nrow(a)) {
    rdmat[j, ] <- scol.matrix.order(a, a[j, ])
  }
  if (direction == 2) {
    rdmat.diff <- rdmat - t(rdmat)
    sel <- which(rdmat.diff > 0)
    rdmat[sel] <- t(rdmat)[sel]
  }
  return(rdmat)
}

gene.specific.p <- function(null.distr, new.d) {
  for (i in 1:length(new.d)) {
    new.d[i] <- pnorm(new.d[i], mean = null.distr[i, 
                                                  1], sd = null.distr[i, 2], lower.tail = TRUE)
  }
  return(new.d)
}

gene.specific.q <- function(new.d) {
  for (i in 1:length(new.d)) {
    new.d[i] <- qnorm(new.d[i], mean = 0, sd = 1, lower.tail = TRUE)
  }
  return(new.d)
}

normscore.row <- function(a) {
  b <- t(apply(a, 1, normal_trafo))
  return(b)
}

#--------------------------------------------------------------------------------
#Actual code
library("fdrtool")
#Sets the directory

args = commandArgs(trailingOnly = TRUE)
infile = args[1]
id = args[2]

mat <- read.table(infile, header = FALSE, sep = '\t')
input = data.matrix(mat)

#Function inputs
use.normal.approx = FALSE
normalization = "standardize"
gene.fdr.plot = FALSE


n <- nrow(input)
m <- ncol(input)
if (use.normal.approx & normalization == "none") 
  normalization <- "standardize"
if (normalization == "normal_score") {
  array <- normscore.row(input)
} else {
  if (normalization == "standardize") {
    array <- normrow(input)
  }
  else {
    if (normalization == "none") {
      array <- input
    }
    else {
      return("wrong normalization method.")
    }
  }
}
orig.array <- array
if (!use.normal.approx) {
  null.distr <- gene.specific.null(array)
} else {
  m.emp <- 2/sqrt(pi) * (m - 1)
  s.emp <- sqrt(m - 1) * 4 * (2 * pi + 3 * sqrt(3) - 9)/3/pi
  null.distr <- cbind(rep(m.emp, nrow(array)), rep(s.emp, 
                                                   nrow(array)))
}
sim.mat <- scol.matrix(array, direction = 1)
d.mat <- sim.mat
d.tmp.mat <- d.mat
for (i in 1:nrow(sim.mat)) d.tmp.mat[i, ] <- gene.specific.p(null.distr, 
                                                             sim.mat[i, ])
for (i in 1:nrow(sim.mat)) d.mat[i, ] <- gene.specific.q(d.tmp.mat[i, 
])
diag(d.mat) <- 0
gene.rel.mat <- matrix(0, nrow = n, ncol = n)
gene.fdr.mat <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  sim.vec <- d.mat[i, ]
  suppressWarnings(t.locfdr <- fdrtool(sim.vec, statistic = "normal", 
                                       plot = gene.fdr.plot, color.figure = TRUE, verbose = FALSE, 
                                       cutoff.method = "pct0", pct0 = 0.75))
  t.row.lfdr <- as.vector(t.locfdr$lfdr)
  gene.fdr.mat[i, ] <- t.row.lfdr
}

#-----------------------------------------------------------------------------
#Saves it as a .csv file
library("MASS")
mat <- (1.0-gene.fdr.mat)
for (i in 1:n){
  for (j in i:n) {
    c_1 = mat[i,j]
    c_2 = mat[j,i]
    mat[i,j] = max(c_1,c_2)
    mat[j,i] = max(c_1,c_2)    
  }
}
file_name <- paste("nlnet_scores_", id, ".csv", sep="")
gene.fdr.mat
write.matrix(mat,file=file_name)

#Now write in what the matrix predicts
place_holder <- mat
for (i in 1:n){
  for (j in i:n) {
    c_1 = mat[i,j]
    c_2 = mat[j,i]
    if ((c_1 <= 0.95 && c_1 >= 0.8) || (c_2 <= 0.95 && c_2 >= 0.8)) {
        place_holder[i,j] = 1
        place_holder[i,j] = 1
    }
    else {
        place_holder[i,j] = 0
        place_holder[i,j] = 0        
    }
  }
}
file_name <- paste("nlnet_", id, ".csv", sep="")
write.matrix(place_holder,file=file_name)
