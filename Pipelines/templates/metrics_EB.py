#This script will take in the matrix from each of the algorithms and the original matrix.
    #It will determine how good the program can detect the values

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import sys

#Firstly let's read in the .csv file
no_prior = np.genfromtxt(sys.argv[1], delimiter = ",")
full_prior = np.genfromtxt(sys.argv[2], delimiter = ",")
genie_prior = np.genfromtxt(sys.argv[3], delimiter = ",")


#Now the actual matrix
original = np.genfromtxt(sys.argv[4], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))

num_genes = sys.argv[5]

#Now converying the matricies into a long list where AUROC comparisons can be made
original = [x for row in original for x in row]
no_prior = [x for row in no_prior for x in row]
full_prior = [x for row in full_prior for x in row]
genie_prior = [x for row in genie_prior for x in row]

#Getting Matthews correlation:
no_prior_matt = metrics.matthews_corrcoef(original,no_prior)
full_prior_matt = metrics.matthews_corrcoef(original,full_prior)
genie_prior_matt = metrics.matthews_corrcoef(original,genie_prior)

#Get a F1 score
no_prior_f1 = metrics.f1_score(original,no_prior)
full_prior_f1 = metrics.f1_score(original,full_prior)
genie_prior_f1 = metrics.f1_score(original,genie_prior)

matthew_name = 'matthew_' + str(num_genes) + '.txt'

with open(matthew_name, 'w') as f:
    f.write('Matthew\'\ correlation:\n')
    f.write('No priors: ' + str(no_prior_matt) + '\n')
    f.write('Full priors: ' + str(full_prior_matt) + '\n')
    f.write('GENIE3 priors: ' + str(genie_prior_matt) + '\n')    
    f.write('\n')
    f.write('F1\'\ score:\n')
    f.write('No priors: ' + str(no_prior_f1) + '\n')
    f.write('Full priors: ' + str(full_prior_f1) + '\n')
    f.write('GENIE3 priors: ' + str(genie_prior_f1) + '\n')   
    f.write('\n')   

#Now calculate the Area under the curve
#Note for now this only calculates the NLNET AUROC, but other metrics can be added once this has been completed. 
fpr_no_prior, tpr_no_prior, thresholds_no_prior = metrics.roc_curve(original,no_prior)
fpr_full_prior, tpr_full_prior, thresholds_full_prior = metrics.roc_curve(original,full_prior)
fpr_genie_prior, tpr_genie_prior, thresholds_genie_prior = metrics.roc_curve(original,genie_prior)

no_prior_score = metrics.auc(fpr_no_prior,tpr_no_prior)
full_prior_score = metrics.auc(fpr_full_prior,tpr_full_prior)
genie_prior_score = metrics.auc(fpr_genie_prior,tpr_genie_prior)


#Plotting the graph
plot_name = plot_name = 'ROCplot_' + str(num_genes) + '.pdf'
plt.plot(fpr_no_prior,tpr_no_prior, color ='r' ,label = 'No prior = %0.3f' %no_prior_score)                #No prior plot
plt.plot(fpr_full_prior,tpr_full_prior, color ='g' ,label = 'Full prior = %0.3f' %full_prior_score)        #Full prior plot
plt.plot(fpr_genie_prior,tpr_genie_prior, color ='b' ,label = 'genie3 prior = %0.3f' %genie_prior_score)   #Genie prior plot
plt.title('ROC curve')                                                                                     #Title
plt.xlabel('FPR')                                                                                          #x label (false positive rate)
plt.ylabel('TPR')                                                                                          #y label (true positive rate)
plt.legend(loc = 'lower right')                                                                            #Legend
plt.savefig(plot_name, format = "pdf")                                                                     #Saves the plot
plt.show()                                                                                                 #Shows the plot

