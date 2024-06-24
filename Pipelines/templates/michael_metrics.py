#This file here will create two sets of metrics:

#1: Get the AUPR curve for each gene set
#2: Get a table showing what genes are shown and what genes are not shown.

#Library:
import sklearn.metrics as metrics
import numpy as np
import sys

#Step 1: Get AUPR curve
#We import in the statistics we need:
pi0_1 = np.genfromtxt(sys.argv[1], delimiter = ",")
original = np.genfromtxt(sys.argv[2], delimiter = ",")
gene_list = open(sys.argv[3], 'r')

#Put it into a workable form
pi0_1_list = [x for row in pi0_1 for x in row]
original_list = [x for row in original for x in row]
gene_list = [gene.rstrip('\n') for gene in gene_list]

#Calculate the AUPR
precision, recall, thresholds = metrics.precision_recall_curve(original_list,pi0_1_list)
metrics.auc(recall,precision)

#Step 2: Now we need a table showing all the gene connections:
#Three lists:
agreed = [] #Connections both found
sim_found = [] #Connections only found by EB
orig_found = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(pi0_1.shape[0]):
    for j in range(i,pi0_1.shape[1]):
        if pi0_1[i][j] == 1:
            sim_found.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed = set(sim_found).intersection(set(orig_found))
agreed = list(agreed)
sim_found = [x for x in sim_found if x not in agreed]
orig_found = [x for x in orig_found if x not in agreed]

#Step 3: Now we put this data into a text file
#Here are all the variables we will use for naming the file
w0 = sys.argv[4]
type = sys.argv[5]
# prior = sys.argv[6]
file_name = 'metrics_' + str(type) + '_' + str(w0) + '.txt'
with open(file_name, 'w') as f:
    f.write('AUPR score: ' + str(metrics.auc(recall,precision)) + '\n\n')
    f.write('Agreed upon gene connections (' + str(len(agreed)) + '):  \n')
    for connection in agreed:
        f.write(connection + '\n')
    f.write('\n')
    f.write('Missed gene connections (' + str(len(orig_found)) + '):  \n')
    for connection in orig_found:
        f.write(connection + '\n')
    f.write('\n')
    f.write('Incorrect gene connections (' + str(len(sim_found)) + '):  \n')
    for connection in sim_found:
        f.write(connection + '\n')
    f.write('\n')


