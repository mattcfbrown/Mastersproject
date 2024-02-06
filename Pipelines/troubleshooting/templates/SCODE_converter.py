#This takes the output of SCODE inference and works out the network connections:

import numpy as np

file = np.loadtxt('/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Data/TestSCODE.txt', delimiter='\t')

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
        