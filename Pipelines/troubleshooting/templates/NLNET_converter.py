#This program will try to extract the information generated from the network inference program

import numpy as np
import itertools
import sys

file = open(sys.argv[1], 'r')

num_genes = sys.argv[2]
num_genes = int(num_genes)

data = []
#Extracts the raw matrix file
for Lines in file.readlines():
    if Lines[0] == '[' or Lines[1] == '[':
        Lines = Lines.replace(' ', '')
        Lines = Lines.replace('[', '')
        Lines = Lines.replace(']', '')
        Lines = Lines.replace('\n', '')
        Lines = Lines.split(",")[1]                
        data.append(Lines)


split = len(data)/num_genes

#Now create the array
rows = []
cur_row = []
for i in range(int(len(data)/split)):
    for j in range(int(split)):
        cur_row.extend(data[i + num_genes*j])
    rows.append(cur_row)
    cur_row = []

a_matrix = np.array(rows)
np.savetxt("matrix_nlnet.csv", a_matrix.astype(int), delimiter=",")
