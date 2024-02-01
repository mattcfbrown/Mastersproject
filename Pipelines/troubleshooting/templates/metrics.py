#This script will take in the matrix from each of the algorithms and the original matrix.
    #It will determine how good the program can detect the values

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import sys

#Firstly let's read in the .csv file
nlnet_Data = np.genfromtxt(sys.argv[1], delimiter = ",")
ni_Data = np.genfromtxt(sys.argv[2], delimiter = ",")
genie_Data = np.genfromtxt(sys.argv[3], delimiter = ",")

#Now the actual matrix
original = np.genfromtxt(sys.argv[4], delimiter=",")
original[0][0] = 0

#Now converying the matricies into a long list where AUROC comparisons can be made
original = [x for row in original for x in row]
nlnet_Data = [x for row in nlnet_Data for x in row]
ni_Data = [x for row in ni_Data for x in row]
genie_Data = [x for row in genie_Data for x in row]

#Getting Matthews correlation:
nlnet_matt = metrics.matthews_corrcoef(original,nlnet_Data)
ni_matt = metrics.matthews_corrcoef(original,ni_Data)
genie_matt = metrics.matthews_corrcoef(original,genie_Data)

with open('matthew.txt', 'w') as f:
    f.write('Matthew\'\ correlation:\n')
    f.write('NLNET: ' + str(nlnet_matt) + '\n')
    f.write('Information measures: ' + str(ni_matt) + '\n')    
    f.write('GENIE3: ' + str(genie_matt) + '\n')

#Now calculate the Area under the curve
#Note for now this only calculates the NLNET AUROC, but other metrics can be added once this has been completed. 
fpr_nlnet, tpr_nlnet, thresholds_nlnet = metrics.roc_curve(original,nlnet_Data)
fpr_ni, tpr_ni, thresholds_ni = metrics.roc_curve(original,ni_Data)
fpr_genie, tpr_genie, thresholds_genie = metrics.roc_curve(original,genie_Data)

nlnet_score = metrics.auc(fpr_nlnet,tpr_nlnet)
ni_score = metrics.auc(fpr_ni,tpr_ni)
genie_score = metrics.auc(fpr_genie,tpr_genie)


#Plotting the graph
plt.plot(fpr_nlnet,tpr_nlnet, color ='r' ,label = 'nlnet = %0.3f' %nlnet_score)                #nlnet plot
plt.plot(fpr_ni,tpr_ni, color ='g' ,label = 'Information_measures = %0.3f' %ni_score)          #ni plot
plt.plot(fpr_genie,tpr_genie, color ='g' ,label = 'GENIE 3 = %0.3f' %ni_score)                 #genie plot
plt.title('ROC curve')                                                                         #Title
plt.xlabel('FPR')                                                                              #x label (false positive rate)
plt.ylabel('TPR')                                                                              #y label (true positive rate)
plt.legend(loc = 'lower right')                                                                #Legend
plt.savefig("ROCplot.pdf", format = "pdf")                                                     #Saves the plot
plt.show()                                                                                     #Shows the plot

