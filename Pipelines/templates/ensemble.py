#This here will take in all the predictions thresholds and convert it into an ensemble:

import numpy as np
import sys

#Import the original matrix (we will use this to determine size)
original = np.genfromtxt(sys.argv[2], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))

num_datasets = len(sys.argv) - 3

#We do this for all the inputs
value = np.zeros((num_genes,num_genes))
for i in range(3,len(sys.argv)):
    data = np.genfromtxt(sys.argv[i], delimiter=",")
    for j in range(num_genes):
        for k in range(num_genes):
            value[j][k] = value[j][k] + (data[j][k]/num_datasets)


id = sys.argv[1]
name = "ensemble_" + id + ".csv"
np.savetxt(name, value, delimiter=",")

