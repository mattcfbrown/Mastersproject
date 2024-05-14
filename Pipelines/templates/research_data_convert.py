#This file takes in the research data which I have taken and converts it into something which can be used

import numpy as np
import sys

data = sys.argv[1]
truth = sys.argv[2]
number = sys.argv[3]
multi = sys.argv[4]

#Generated from:
#https://stackoverflow.com/questions/71758470/numpy-reading-a-csv-file-with-string-header-and-row-names
Get_names = np.genfromtxt(data, delimiter= ',', skip_header=0,dtype=None)
#Gets the gene names
gene_names = [name.decode('UTF-8') for name in Get_names[:,0]][1:]

#Puts the read data into a numpy array, which we can return as a .csv file
test = np.genfromtxt(data, delimiter= ',', skip_header=1)
test = test[:,1:]
np.savetxt("expression" + str(number) + "_" + str(multi) + ".txt", test.astype(float), delimiter="\t")



ground_truth = np.genfromtxt(truth, delimiter= ',', skip_header = 1, dtype=str)
#We now generate a .csv file which can be used as the ground truth network
num_genes = len(gene_names)
truth_matrix = np.zeros((num_genes,num_genes))
for i in range(len(ground_truth)):
    gene_1 = gene_names.index(ground_truth[i,0])
    gene_2 = gene_names.index(ground_truth[i,1])
    truth_matrix[gene_1][gene_2] = 1
    truth_matrix[gene_2][gene_1] = 1

np.savetxt("ground_truth" + str(number) + "_" + str(multi) + ".csv", truth_matrix.astype(int), delimiter=",")


#We also want a list of genes
file = open("genes" + str(number) + "_" + str(multi) + ".txt", 'w')
for genes in gene_names:
    file.write(genes + "\n")
file.close()
