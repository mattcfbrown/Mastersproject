#This program will try to extract the information generated from the network inference program

import numpy as np

file = open('/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/troubleshooting/results/Information_Measures/outfile_NI.txt', 'r')

#Extracts the raw matrix file
for Lines in file.readlines():
    if 'Bool' in Lines:
        matrix = Lines


#Now let's create a matrix from this data
#Removes data not needed
matrix = matrix.replace('Bool', '')
matrix = matrix.replace('[', '')
matrix = matrix.replace(']', '')
matrix = matrix.replace(' ', '')
matrix = matrix.replace('\n', '')


#Now create the array
rows = []
cur_row = []
for i in range(len(matrix)):
    if matrix[i] == ';':
        rows.append(cur_row)
        cur_row = []
    else:
        cur_row.append(int(matrix[i]))

a_matrix = np.array(rows)
print(a_matrix)
