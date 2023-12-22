import matplotlib.pyplot as plt
import numpy as np

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

genes = 25

#This is for information measures
with open("/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Julia_programs/outfile_NI.txt", 'r+') as f:
    for line in f:
        if line.startswith("Bool"):
            IM = line
            break

#IM = IM[4:]
#IM = IM.replace("[", "")
#IM = IM.replace("]", "")
#IM = np.matrix(IM)

#Now for NLNET
NLNET = []
found = 0
with open("/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Julia_programs/outfile_NLNET.txt", 'r+') as f:
    for line in f:
        if line.startswith("$graph"):
            found = 1
        if found == 1 and line != '\n':
            NLNET.append(line)
        if found == 1 and line == '\n':
            break


NLNET_matrix = np.zeros(shape = (genes,genes))
#Now clean it up a bit
for line in NLNET:
    #Gets all the data which is needed
    if '[' in line:
        for i in range(len(line)):
            if line[i:i+2] == '--':
                if line[i-2] != '' and i >2:
                    back = int(line[i-2:i])
                if line[i+3] != '':
                    front = int(line[i+2:i+4])
                NLNET_matrix[back,front] = 1
                NLNET_matrix[front,back] = 1

#This converts it to a useable form we can work with
NLNET_list = list(map(int,NLNET_matrix.flatten()))




# generate 2 class dataset
#y = [1,1,1,1,1,0,0,0,0,0]
#X = [1,1,1,0,0,0,0,1,1,0]

#auc = roc_auc_score(y, X)
#print('AUC: %.3f' % auc)

#ns_fpr, ns_tpr, _ = roc_curve(y, X)

#plt.plot(ns_fpr, ns_tpr, linestyle='--', label='No Skill')
# axis labels
#plt.xlabel('False Positive Rate')
#plt.ylabel('True Positive Rate')
# show the legend
#plt.legend()
# show the plot
#plt.show()