#This script will take in the matrix from each of the algorithms and the original matrix.
    #It will determine how good the program can detect the values

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import sys

#Firstly let's read in the .csv file
nlnet_Data = np.genfromtxt(sys.argv[1], delimiter = ",")

#Now the actual matrix
original = np.genfromtxt(sys.argv[2], delimiter=",")
original[0][0] = 0

#Now converying the matricies into a long list where AUROC comparisons can be made
original = [x for row in original for x in row]
nlnet_Data = [x for row in nlnet_Data for x in row]

#nlnet_Data = [1 if x>1 else x for x in nlnet_Data]

#Now calculate the Area under the curve
#Note for now this only calculates the NLNET AUROC, but other metrics can be added once this has been completed. 
fpr, tpr, thresholds = metrics.roc_curve(original,nlnet_Data)
nlnet_score = metrics.auc(fpr,tpr)

#Plotting the graph
plt.plot(fpr,tpr, label = 'nlnet = %0.3f' %nlnet_score)                #The plot itself
plt.title('ROC curve')                                                 #Title
plt.xlabel('FPR')                                                      #x label (false positive rate)
plt.ylabel('TPR')                                                      #y label (true positive rate)
plt.legend(loc = 'lower right')                                        #Legend
plt.savefig("ROCplot.pdf", format = "pdf")                             #Saves the plot
plt.show()                                                             #Shows the plot

