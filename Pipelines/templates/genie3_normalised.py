#This script converts the output of GENIE3 into a csv file

import numpy as np
import sys
import re
import math

infile = sys.argv[1]
num_genes = sys.argv[2]
num_genes = int(num_genes)

#This produces a list containing all the values
data = []
with open(infile, 'r') as fp:
    lines = fp.readlines()

for line in lines[1:]:
    #Gets the location of the genes
    location =[loc.start() for loc in re.finditer('Gene', line)]
    #This gives us the gene number
    location = [x+4 for x in location]
    location.append(0)
    vals = []
    for loc in location:
        current = loc
        gene_ID = ''
        while line[current] != ' ' and current < len(line)-1:
            gene_ID = gene_ID + line[current]
            current = current + 1
        location[2] = location[1] + len(gene_ID) + 1   
        vals.append(gene_ID)
    data.append(vals)


matrix = np.zeros((num_genes,num_genes))

max = float(data[0][2])
min = float(data[math.perm(num_genes,2)-1][2])
print(max)
print(min)


#now puts into matrix
for x in data:
    val1 = int(x[0]) - 1
    val2 = int(x[1]) - 1
    matrix[val1][val2] = (float(x[2])-min)/(max-min)


np.savetxt("matrix_genie_normalised.csv", matrix, delimiter=",")
