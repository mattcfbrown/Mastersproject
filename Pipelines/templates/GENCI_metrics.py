#Here I will get the metrics for my files

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import sys
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join

#We firstly read in all the inference files which we have
path = sys.argv[1]
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#We then read in both EB files
eb_full = sys.argv[2]
eb_ten = sys.argv[3]
onlyfiles.append(eb_full)
onlyfiles.append(eb_ten)

#We now read in the ground truth network
original = sys.argv[4]
original = np.genfromtxt(original,delimiter=',')
num_genes = original.shape[0]
original[0][0] = 0

#This is the ID value
id = sys.argv[5]

#Ensemble reading in

ensemble_full = sys.argv[6]
ensemble_ten = sys.argv[7]
onlyfiles.append(ensemble_full)
onlyfiles.append(ensemble_ten)

arrays = []    #Holds all the arrays used
names = []     #Holds the names of the methods
AUPR = []      #Holds AUPR scores
found = []     #Holds all found
missed = []    #Holds the missed
incorrect = [] #Holds the incorrect
num_connections = num_genes*num_genes - num_genes

#Making a quick gene_list
gene_list = []
for i in range(num_genes):
    gene_list.append("G" + str(i))

#For each file we wish to do the following:
for method in onlyfiles:
    #Firstly get the name of the method
    name = method.replace("GRN_", "")
    name = name.replace(".csv", "")
    name = name.replace("_" + id, "")
    name = name.replace("genci_input_folder/all/final_list", "ensemble_all")
    name = name.replace("genci_input_folder/top_10/final_list", "ensemble_ten")
    print(name)
    names.append(name)

    #Secondly, we want to create a file 
    if name == 'EB_full':
        matrix = np.genfromtxt(eb_full, delimiter=',',
                        dtype=None)
    elif name == 'EB_ten':
        matrix = np.genfromtxt(eb_ten, delimiter=',',
                        dtype=None)
    else:
        if name == 'ensemble_all':
            file = np.genfromtxt(ensemble_full, delimiter=',',
                        dtype=None)
        elif name == 'ensemble_ten':
            file = np.genfromtxt(ensemble_ten, delimiter=',',
                        dtype=None)
        else:
            file = np.genfromtxt(path + '/' + method, delimiter=',',
                        dtype=None)
        #We only want the top 10% of connections
        matrix = np.zeros(shape=(num_genes,num_genes))
        hit = 0
        for i in range(len(file)):
            gene1 = int(file[i][0][1:])
            gene2 = int(file[i][1][1:])
            #Ad Hoc solution to an error I am getting
            if num_genes == 20:
                gene1 = gene1-1
                gene2 = gene2-1
            if matrix[gene1][gene2] == 0:
                matrix[gene1][gene2] = 1
                matrix[gene2][gene1] = 1
                hit = hit + 1
            if hit == round(num_connections/20):
                break
    
    #Now get the AUPR scores
    matrix_list = [x for row in matrix for x in row]
    original_list = [x for row in original for x in row]
    precision, recall, thresholds = metrics.precision_recall_curve(original_list,matrix_list)
    AUPR.append(metrics.auc(recall,precision))

    #Now get the information about the missing and found gene connections
    found_current = [] #Connections found
    missed_current = [] #Connections missed
    incorrect_current = [] #Connections incorrectly guessed

    for i in range(matrix.shape[0]):
        for j in range(i,matrix.shape[1]):
            if matrix[i][j] == 1 and original[i][j] == 1:
                found_current.append(gene_list[i] + ' --- ' + gene_list[j])
            if matrix[i][j] == 1 and original[i][j] == 0:
                incorrect_current.append(gene_list[i] + ' --- ' + gene_list[j]) 
            if matrix[i][j] == 0 and original[i][j] == 1:
                missed_current.append(gene_list[i] + ' --- ' + gene_list[j])

    found.append(found_current)
    incorrect.append(incorrect_current)
    missed.append(missed_current)

#We now read all which is above into a file
file_name = 'metrics_' + str(id) + '.txt'
with open(file_name, 'w') as f:
    #Firstly write in AUPR scores
    f.write('AUPR score: ' + '\n\n')
    for i in range(len(names)):
        f.write(names[i] + ': ')
        f.write(str(AUPR[i]) + '\n')
    f.write('\n')
    #Now time to write in the missed
    for i in range(len(names)):
        f.write(names[i] + ' connection information\n\n')
        f.write('FOUND:\n')
        for connection in found[i]:
            f.write(connection)
            f.write('\n')
        f.write('\n')
        f.write('MISSED:\n')
        for connection in missed[i]:
            f.write(connection)
            f.write('\n')
        f.write('\n')
        f.write('INCORRECT:\n')
        for connection in incorrect[i]:
            f.write(connection)
            f.write('\n')
        f.write('\n\n')
        


