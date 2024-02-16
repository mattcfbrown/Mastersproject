#This takes the output of SCODE inference and works out the network connections:

from math import floor
import numpy as np
import sys

input = sys.argv[1]
file = np.loadtxt(input, delimiter='\t')

#What we are going to do here is rank the absolute values produced.
#2D ARRAY:
#COL 1 = Gene 1
#Col 2 = Gene 2
#Col 3 = score

num_genes = file.shape[0]
scores = np.zeros(shape=(num_genes*num_genes, 3))

current = 0
for i in range(num_genes):
    for j in range(num_genes):
        scores[current][0] = i
        scores[current][1] = j
        scores[current][2] = abs(file[i][j])
        current = current + 1


#We now sort the scores array to find the highest absolute values
scores = scores[scores[:, 2].argsort()[::-1]]

#We create a matrix, to which we can store our predictions
matrix = np.zeros(shape = (num_genes,num_genes))
#The total number of possible gene to gene connections which can exist
total_con = num_genes*num_genes
i = 0

while i <= int(floor(total_con/2)-1):
    gene1 = int(scores[i][0])
    gene2 = int(scores[i][1])
    matrix[gene1][gene2] = 1
    i = i + 1

np.savetxt("matrix_SCODE.csv", matrix, delimiter=",")

        