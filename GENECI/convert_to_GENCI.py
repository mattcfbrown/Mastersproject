#This aims to try and convert some of the data I already have into a compatable format:

import numpy as np
import pandas as pd

input = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Ensemble_testing/Data/25genes_250cells.txt"
data = np.genfromtxt(input, delimiter="\t")

gene_names = []
#Gene names
for i in range(data.shape[0]):
    name = "G" + str(i)
    gene_names.append(name)

cell_names = []
#Cell names
for j in range(data.shape[1]):
    name = "C" + str(j)
    cell_names.append(name)

df = pd.DataFrame(data, index=gene_names, columns=cell_names)
print(df)
df.to_csv('test.csv', index=True, header=True, sep=',')
