#This script aims to do two things
#Conver GENIE3 output into a csv file
#Convert the GENIE3 into a ground truth network (i.e see what GENIE3 would predict)

import numpy as np
import sys
import re
import math

infile = sys.argv[1]
type = sys.argv[2]

#This produces a list containing all the values
data = []
with open(infile, 'r') as fp:
    lines = fp.readlines()

# num_genes = (1 + ( float(1 + 4*len(lines[1:])) )^0.5 )/2
num_genes = 25
num_genes = int(num_genes)

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

max = float(data[0][2])
min = float(data[math.perm(num_genes,2)-1][2])

matrix = np.zeros((num_genes,num_genes))

#now puts into matrix
for x in data:
    val1 = int(x[0]) - 1
    val2 = int(x[1]) - 1
    matrix[val1][val2] = (float(x[2]) - min)/(max-min)

#This here is simply the Matrix which will be used as priors
np.savetxt("matrix_genie_eb_" + type + ".csv", matrix, delimiter=",")

#We now will make the GENIE predicted matrix
genie_matrix = np.zeros((num_genes,num_genes))
threshold = 0.1
for x in data:
    if float(x[2]) < threshold:
        break
    val1 = int(x[0])-1
    val2 = int(x[1])-1
    genie_matrix[val1][val2] = 1
    genie_matrix[val2][val1] = 1

np.savetxt("genie_predict_" + type + ".csv", genie_matrix, delimiter=",")