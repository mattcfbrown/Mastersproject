#This script converts the output of GENIE3 into a csv file

import numpy as np
import sys
import re

file = open('/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/troubleshooting/results/genie3/genie.txt','r')
num_genes = 9

#This produces a list containing all the values
data = []
for Lines in file.readlines():
    if 'Gene' in Lines and 'regulatoryGene' not in Lines:
        #Gets the location of the genes
        location =[loc.start() for loc in re.finditer('Gene', Lines)]
        #This gives us the gene number
        location = [x+4 for x in location]
        vals = []
        for loc in location:
            current = loc
            gene_ID = ''
            while Lines[current] != ' ':
                gene_ID = gene_ID + Lines[current]
                current = current + 1
            vals.append(gene_ID)
        data.append(vals)

matrix = np.zeros((num_genes,num_genes))

#now puts into matrix
for x in data:
    val1 = int(x[0]) - 1
    val2 = int(x[1]) - 1
    matrix[val1][val2] = 1
    matrix[val2][val1] = 1

np.savetxt("matrix_genie.csv", matrix, delimiter=",")
