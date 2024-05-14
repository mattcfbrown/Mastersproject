#This script will take in the threshold values and output some metrics

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import sys

#Firstly let's read in the .csv file
first_threshold = np.genfromtxt(sys.argv[1], delimiter = ",")
second_threshold = np.genfromtxt(sys.argv[2], delimiter = ",")
third_threshold = np.genfromtxt(sys.argv[3], delimiter = ",")
fourth_threshold = np.genfromtxt(sys.argv[4], delimiter = ",")
fifth_threshold = np.genfromtxt(sys.argv[5], delimiter = ",")
#Now the original
original = np.genfromtxt(sys.argv[6], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))

threshold = sys.argv[7].split(',')
type = sys.argv[8]
method = sys.argv[9]
num_cells = sys.argv[10]

#Now converying the matricies into a long list where AUROC comparisons can be made
original = [x for row in original for x in row]
first_threshold = [x for row in first_threshold for x in row]
second_threshold = [x for row in second_threshold for x in row]
third_threshold = [x for row in third_threshold for x in row]
fourth_threshold = [x for row in fourth_threshold for x in row]
fifth_threshold = [x for row in fifth_threshold for x in row]

#Getting Matthews correlation:
first_matt = metrics.matthews_corrcoef(original,first_threshold)
second_matt = metrics.matthews_corrcoef(original,second_threshold)
third_matt = metrics.matthews_corrcoef(original,third_threshold)
fourth_matt = metrics.matthews_corrcoef(original,fourth_threshold)
fifth_matt = metrics.matthews_corrcoef(original,fifth_threshold)

#Get a F1 score
first_f1 = metrics.f1_score(original,first_threshold)
second_f1 = metrics.f1_score(original,second_threshold)
third_f1 = metrics.f1_score(original,third_threshold)
fourth_f1 = metrics.f1_score(original,fourth_threshold)
fifth_f1 = metrics.f1_score(original,fifth_threshold)

matthew_name = 'matthew_' + str(num_genes) + '_' + str(num_cells) + '_' + type + '_'  + method +'.txt'
with open(matthew_name, 'w') as f:
    f.write('Matthew\'\ correlation:\n')
    f.write(str(threshold[0]) + ' :' + str(first_matt) + '\n')
    f.write(str(threshold[1]) + ' :' + str(second_matt) + '\n')
    f.write(str(threshold[2])+ ' :' + str(third_matt) + '\n')
    f.write(str(threshold[3]) + ' :' + str(fourth_matt) + '\n')    
    f.write(str(threshold[4]) + ' :' + str(fifth_matt) + '\n')   
    f.write('\n')
    f.write('F1\'\ score:\n')
    f.write(str(threshold[0]) + ' :' + str(first_f1) + '\n')
    f.write(str(threshold[1])+ ' :' + str(second_f1) + '\n')
    f.write(str(threshold[2]) + ' :' + str(third_f1) + '\n')
    f.write(str(threshold[3]) + ' :' + str(fourth_f1) + '\n')    
    f.write(str(threshold[4]) + ' :' + str(fifth_f1) + '\n')      
    f.write('\n')   

#Now calculate the Area under the curve
#Note for now this only calculates the NLNET AUROC, but other metrics can be added once this has been completed. 
precision_first, recall_first, thresholds_first = metrics.precision_recall_curve(original,first_threshold)
precision_second, recall_second, thresholds_second = metrics.precision_recall_curve(original,second_threshold)
precision_third, recall_third, thresholds_third = metrics.precision_recall_curve(original,third_threshold)
precision_fourth, recall_fourth, thresholds_fourth = metrics.precision_recall_curve(original,fourth_threshold)
precision_fifth, recall_fifth, thresholds_fifth = metrics.precision_recall_curve(original,fifth_threshold)

first_score = metrics.auc(recall_first,precision_first)
second_score = metrics.auc(recall_second,precision_second)
third_score = metrics.auc(recall_third,precision_third)
fourth_score = metrics.auc(recall_fourth,precision_fourth)
fifth_score = metrics.auc(recall_fifth,precision_fifth)

#Plotting the graph
plot_name = plot_name = 'ROCplot_' + str(num_genes) + '_' + str(num_cells) + '_' + str(type) + '_'  + str(method) + '.pdf'
plt.plot(precision_first,recall_first, color ='r' ,label = str(threshold[0]) + ' = %0.3f' %first_score)                
plt.plot(precision_second,recall_second, color ='g' ,label = str(threshold[1]) + ' = %0.3f' %second_score)        
plt.plot(precision_third,recall_third, color ='b' ,label = str(threshold[2]) + ' = %0.3f' %third_score)   
plt.plot(precision_fourth,recall_fourth, color ='c' ,label = str(threshold[3]) + ' = %0.3f' %fourth_score)    
plt.plot(precision_fifth,recall_fifth, color ='m' ,label = str(threshold[4]) + ' = %0.3f' %fifth_score)   
plt.title('ROC curve')                                                                                     #Title
plt.xlabel('precision')                                                                                    #x label (false positive rate)
plt.ylabel('recall')                                                                                       #y label (true positive rate)
plt.legend(loc = 'lower right')                                                                            #Legend
plt.savefig(plot_name, format = "pdf")                                                                     #Saves the plot
plt.show()             