#Simple, I just wish to normalise the values

import numpy as np
from math import log
import sys

path = sys.argv[1]
data = np.genfromtxt(path, dtype=None, delimiter='\t')

total = data.sum()
L = (1/(data.size/len(data)))*total #Getting L value
y0 = 1                              #Getting y0 value
columnsum = data.sum(axis=0)
columnsum = columnsum/L
index_locations = np.where(columnsum == 0)[0]
new = np.zeros((len(data),len(columnsum)-len(index_locations)))
#Row then column
#Rows
for i in range(len(data)):
    #Columns
    index = 0
    for j in range(len(columnsum)):
        if j not in index_locations:
            new[i][index] = log((data[i][j]/columnsum[j]) + 1)
            index = index + 1

cells = sys.argv[2]

name = 'normalised_file_' + cells + '.txt'

np.savetxt(name,new,delimiter='\t')


