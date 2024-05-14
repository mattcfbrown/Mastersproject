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
SCODE_Data = np.genfromtxt(sys.argv[4], delimiter = ",")

#Now the actual matrix
original = np.genfromtxt(sys.argv[5], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))

#Now converying the matricies into a long list where AUROC comparisons can be made
original = [x for row in original for x in row]
nlnet_Data = [x for row in nlnet_Data for x in row]
ni_Data = [x for row in ni_Data for x in row]
genie_Data = [x for row in genie_Data for x in row]
SCODE_Data = [x for row in SCODE_Data for x in row]

#Getting Matthews correlation:
nlnet_matt = metrics.matthews_corrcoef(original,nlnet_Data)
ni_matt = metrics.matthews_corrcoef(original,ni_Data)
genie_matt = metrics.matthews_corrcoef(original,genie_Data)
SCODE_matt = metrics.matthews_corrcoef(original,SCODE_Data)

#Get a F1 score
nlnet_f1 = metrics.f1_score(original,nlnet_Data)
ni_f1 = metrics.f1_score(original,ni_Data)
genie_f1 = metrics.f1_score(original,genie_Data)
SCODE_f1 = metrics.f1_score(original,SCODE_Data)

matthew_name = 'matthew_' + str(num_genes) + '.txt'

with open(matthew_name, 'w') as f:
    f.write('Matthew\'\ correlation:\n')
    f.write('NLNET: ' + str(nlnet_matt) + '\n')
    f.write('Information measures: ' + str(ni_matt) + '\n')    
    f.write('GENIE3: ' + str(genie_matt) + '\n')
    f.write('SCODE: ' + str(SCODE_matt) + '\n')
    f.write('\n')
    f.write('F1\'\ score:\n')
    f.write('NLNET: ' + str(nlnet_f1) + '\n')
    f.write('Information measures: ' + str(ni_f1) + '\n')    
    f.write('GENIE3: ' + str(genie_f1) + '\n')
    f.write('SCODE: ' + str(SCODE_f1) + '\n')

#Now calculate the Area under the curve
#Note for now this only calculates the NLNET AUROC, but other metrics can be added once this has been completed. 
precision_nlent, recall_nlnet, thresholds_nlnet = metrics.precision_recall_curve(original,nlnet_Data)
precision_ni, recall_ni, thresholds_ni = metrics.precision_recall_curve(original,ni_Data)
precision_genie, recall_genie, thresholds_genie = metrics.precision_recall_curve(original,genie_Data)
precision_SCODE, recall_SCODE, thresholds_SCODE = metrics.precision_recall_curve(original,SCODE_Data)

nlnet_score = metrics.auc(recall_nlnet,precision_nlent)
ni_score = metrics.auc(recall_ni,precision_ni)
genie_score = metrics.auc(recall_genie,precision_genie)
scode_score = metrics.auc(recall_SCODE,precision_SCODE)


#Plotting the graph
plot_name = 'ROCplot_' + str(num_genes) + '.pdf'
plt.plot(recall_nlnet,precision_nlent, color ='r' ,label = 'nlnet = %0.3f' %nlnet_score)                #nlnet plot
plt.plot(recall_ni,precision_ni, color ='g' ,label = 'Information_measures = %0.3f' %ni_score)          #ni plot
plt.plot(recall_genie,precision_genie, color ='b' ,label = 'GENIE 3 = %0.3f' %genie_score)              #genie plot
plt.plot(recall_SCODE,precision_SCODE, color ='k' ,label = 'SCODE = %0.3f' %scode_score)                #SCODE plot
plt.title('ROC curve')                                                                         #Title
plt.xlabel('FPR')                                                                              #x label (false positive rate)
plt.ylabel('TPR')                                                                              #y label (true positive rate)
plt.legend(loc = 'lower right')                                                                #Legend
plt.savefig(plot_name, format = "pdf")                                                     #Saves the plot
plt.show()                                                                                     #Shows the plot

