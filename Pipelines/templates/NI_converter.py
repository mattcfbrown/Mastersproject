#This program will try to extract the information generated from the network inference program

import numpy as np
import sys

file = open(sys.argv[1], 'r')
num_genes = sys.argv[2]
num_genes = int(num_genes)


found = 0
matrix = ''
#Extracts the raw matrix file
for Lines in file.readlines():
    if 'Bool' in Lines:
        found = 1
    elif found == 1 and ('0' or '1' in Lines):
        matrix = matrix + Lines


#Now let's create a matrix from this data
#Removes data not needed
matrix = matrix.replace(' ', '')
matrix = matrix.replace('\n', '')


#Now create the array
rows = []
cur_row = []
for i in range(len(matrix)):
    if (i+1) % num_genes == 0:
        cur_row.append(int(matrix[i]))
        rows.append(cur_row)
        cur_row = []
    else:
        cur_row.append(int(matrix[i]))


a_matrix = np.array(rows)
np.savetxt("matrix_NI.csv", a_matrix.astype(int), delimiter=",")
