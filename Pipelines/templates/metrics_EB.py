#This script will take in the matrix from each of the algorithms and the original matrix.
    #It will determine how good the program can detect the values

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import sys

#Firstly let's read in the .csv file
no_prior = np.genfromtxt(sys.argv[1], delimiter = ",")
full_prior = np.genfromtxt(sys.argv[2], delimiter = ",")
# genie_prior = np.genfromtxt(sys.argv[3], delimiter = ",")
# nlnet_prior = np.genfromtxt(sys.argv[4], delimiter = ",")
# genie_normal = np.genfromtxt(sys.argv[5], delimiter = ",")
# ni_normal = np.genfromtxt(sys.argv[6], delimiter = ",")

#Now the actual matrix
original = np.genfromtxt(sys.argv[3], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))

num_genes = sys.argv[4]
num_genes = int(num_genes)
num_cells = sys.argv[5]

#Now converying the matricies into a long list where AUROC comparisons can be made
original_list = [x for row in original for x in row]
no_prior_list = [x for row in no_prior for x in row]
full_prior_list = [x for row in full_prior for x in row]
# genie_prior_list = [x for row in genie_prior for x in row]
# nlnet_prior_list = [x for row in nlnet_prior for x in row]
# genie_normal_list = [x for row in genie_normal for x in row]
# ni_normal_list = [x for row in ni_normal for x in row]

#Getting Matthews correlation:
no_prior_matt = metrics.matthews_corrcoef(original_list,no_prior_list)
full_prior_matt = metrics.matthews_corrcoef(original_list,full_prior_list)
# genie_prior_matt = metrics.matthews_corrcoef(original_list,genie_prior_list)
# nlnet_prior_matt = metrics.matthews_corrcoef(original_list,nlnet_prior_list)
# genie_normal_matt = metrics.matthews_corrcoef(original_list,genie_normal_list)
# ni_normal_matt = metrics.matthews_corrcoef(original_list,ni_normal_list)

#Get a F1 score
no_prior_f1 = metrics.f1_score(original_list,no_prior_list)
full_prior_f1 = metrics.f1_score(original_list,full_prior_list)
# genie_prior_f1 = metrics.f1_score(original_list,genie_prior_list)
# nlnet_prior_f1 = metrics.f1_score(original_list,nlnet_prior_list)
# genie_normal_f1 = metrics.f1_score(original_list,genie_normal_list)
# ni_normal_f1 = metrics.f1_score(original_list,ni_normal_list)

matthew_name = 'matthew_' + str(num_genes) + '_' + str(num_cells) + '.txt'
 
#Now calculate the Area under the curve
#Note for now this only calculates the NLNET AUROC, but other metrics can be added once this has been completed. 
precision_no_prior, recall_no_prior, thresholds_no_prior = metrics.precision_recall_curve(original_list,no_prior_list)
precision_full_prior, recall_full_prior, thresholds_full_prior = metrics.precision_recall_curve(original_list,full_prior_list)
# precision_genie, recall_genie, thresholds_genie = metrics.precision_recall_curve(original_list,genie_prior_list)
# precision_nlent, recall_nlnet, thresholds_nlnet = metrics.precision_recall_curve(original_list,nlnet_prior_list)
# precision_genie_normal, recall_genie_normal, thresholds_genie_normal = metrics.precision_recall_curve(original_list,genie_normal_list)
# precision_ni, recall_ni, thresholds_ni = metrics.precision_recall_curve(original_list,ni_normal_list)

no_prior_score = metrics.auc(recall_no_prior,precision_no_prior)
full_prior_score = metrics.auc(recall_full_prior,precision_full_prior)
# genie_prior_score = metrics.auc(recall_genie,precision_genie)
# nlnet_prior_score = metrics.auc(recall_nlnet,precision_nlent)
# genie_normal_score = metrics.auc(recall_genie_normal,precision_genie_normal)
# ni_normal_score = metrics.auc(recall_ni,precision_ni)

#Plotting the graph
plot_name = plot_name = 'ROCplot_' + str(num_genes) + '_' + str(num_cells) + '.pdf'
plt.plot(recall_no_prior,precision_no_prior, color ='r' ,label = 'No prior = %0.3f' %no_prior_score)                 #No prior plot
plt.plot(recall_full_prior,precision_full_prior, color ='g' ,label = 'Full prior = %0.3f' %full_prior_score)         #Full prior plot
# plt.plot(recall_genie,precision_genie, color ='b' ,label = 'genie3 prior = %0.3f' %genie_prior_score)    #Genie prior plot
# plt.plot(recall_nlnet,precision_nlent, color ='c' ,label = 'nlnet prior = %0.3f' %nlnet_prior_score)     #nlnet prior plot
# plt.plot(recall_genie_normal,precision_genie_normal, color ='m' ,label = 'GENIE normal = %0.3f' %genie_normal_score) #GENIE normal plot
# plt.plot(recall_ni,precision_ni, color ='y' ,label = 'NI normal = %0.3f' %ni_normal_score)             #NI normal plot
plt.title('ROC curve')                                                                                      #Title
plt.xlabel('FPR')                                                                                           #x label (false positive rate)
plt.ylabel('TPR')                                                                                           #y label (true positive rate)
plt.legend(loc = 'lower right')                                                                             #Legend
plt.savefig(plot_name, format = "pdf")                                                                      #Saves the plot
plt.show()                                                                                                  #Shows the plot


#We now create a gene list for this data, to work out the connected edges
gene_list = []
#We create our gene list
for i in range(num_genes):
    gene_list.append('T' + str(i))

#We now write code for each situation
#No priors
agreed_no = [] #Connections both found
sim_found_no = [] #Connections only found by EB
orig_found_no = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if no_prior [i][j] == 1:
            sim_found_no.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_no.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_no = set(sim_found_no).intersection(set(orig_found_no))
agreed_no = list(agreed_no)
sim_found_no = [x for x in sim_found_no if x not in agreed_no]
orig_found_no = [x for x in orig_found_no if x not in agreed_no]

#Full priors
agreed_full = [] #Connections both found
sim_found_full = [] #Connections only found by EB
orig_found_full = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if full_prior [i][j] == 1:
            sim_found_full.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_full.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_full = set(sim_found_full).intersection(set(orig_found_full))
agreed_full = list(agreed_full)
sim_found_full = [x for x in sim_found_full if x not in agreed_full]
orig_found_full = [x for x in orig_found_full if x not in agreed_full]

#Genie normal
# agreed_g_norm = [] #Connections both found
# sim_found_g_norm = [] #Connections only found by EB
# orig_found_g_norm = [] #Connections not found by EB

# #Let's get the simulated and original ones
# for i in range(num_genes):
#     for j in range(i,num_genes):
#         if genie_normal [i][j] == 1:
#             sim_found_g_norm.append(gene_list[i] + ' --- ' + gene_list[j])
#         if original[i][j] == 1:
#             orig_found_g_norm.append(gene_list[i] + ' --- ' + gene_list[j])            

# #We get agreement
# agreed_g_norm = set(sim_found_g_norm).intersection(set(orig_found_g_norm))
# agreed_g_norm = list(agreed_g_norm)
# sim_found_g_norm = [x for x in sim_found_g_norm if x not in agreed_g_norm]
# orig_found_g_norm = [x for x in orig_found_g_norm if x not in agreed_g_norm]

# #Nlnet prior
# agreed_nlnet = [] #Connections both found
# sim_found_nlnet = [] #Connections only found by EB
# orig_found_nlnet = [] #Connections not found by EB

# #Let's get the simulated and original ones
# for i in range(num_genes):
#     for j in range(i,num_genes):
#         if nlnet_prior [i][j] == 1:
#             sim_found_nlnet.append(gene_list[i] + ' --- ' + gene_list[j])
#         if original[i][j] == 1:
#             orig_found_nlnet.append(gene_list[i] + ' --- ' + gene_list[j])            

# #We get agreement
# agreed_nlnet = set(sim_found_nlnet).intersection(set(orig_found_nlnet))
# agreed_nlnet = list(agreed_nlnet)
# sim_found_nlnet = [x for x in sim_found_nlnet if x not in agreed_nlnet]
# orig_found_nlnet = [x for x in orig_found_nlnet if x not in agreed_nlnet]

# #Genie prior
# agreed_g_prior = [] #Connections both found
# sim_found_g_prior = [] #Connections only found by EB
# orig_found_g_prior = [] #Connections not found by EB

# #Let's get the simulated and original ones
# for i in range(num_genes):
#     for j in range(i,num_genes):
#         if nlnet_prior [i][j] == 1:
#             sim_found_g_prior.append(gene_list[i] + ' --- ' + gene_list[j])
#         if original[i][j] == 1:
#             orig_found_g_prior.append(gene_list[i] + ' --- ' + gene_list[j])            

# #We get agreement
# agreed_g_prior = set(sim_found_g_prior).intersection(set(orig_found_g_prior))
# agreed_g_prior = list(agreed_g_prior)
# sim_found_g_prior = [x for x in sim_found_g_prior if x not in agreed_g_prior]
# orig_found_g_prior = [x for x in orig_found_g_prior if x not in agreed_g_prior]

# #NI normal
# agreed_ni = [] #Connections both found
# sim_found_ni = [] #Connections only found by EB
# orig_found_ni = [] #Connections not found by EB

# #Let's get the simulated and original ones
# for i in range(num_genes):
#     for j in range(i,num_genes):
#         if ni_normal [i][j] == 1:
#             sim_found_ni.append(gene_list[i] + ' --- ' + gene_list[j])
#         if original[i][j] == 1:
#             orig_found_ni.append(gene_list[i] + ' --- ' + gene_list[j])            

# #We get agreement
# agreed_ni = set(sim_found_ni).intersection(set(orig_found_ni))
# agreed_ni = list(agreed_ni)
# sim_found_ni = [x for x in sim_found_ni if x not in agreed_ni]
# orig_found_ni = [x for x in orig_found_ni if x not in agreed_ni]


#Now that we have our 6, we can now write it into text files
with open(matthew_name, 'w') as f:
    f.write('Matthew\'\ correlation:\n')
    f.write('No priors: ' + str(no_prior_matt) + '\n')
    f.write('Full priors: ' + str(full_prior_matt) + '\n')
    # f.write('GENIE3 priors: ' + str(genie_prior_matt) + '\n')
    # f.write('NLNET priors: ' + str(nlnet_prior_matt) + '\n')    
    # f.write('GENIE normal: ' + str(genie_normal_matt) + '\n')
    # f.write('NI normal: ' + str(ni_normal_matt) + '\n')      
    f.write('\n')
    f.write('F1\'\ score:\n')
    f.write('No priors: ' + str(no_prior_f1) + '\n')
    f.write('Full priors: ' + str(full_prior_f1) + '\n')
    # f.write('GENIE3 priors: ' + str(genie_prior_f1) + '\n')
    # f.write('nlnet priors: ' + str(nlnet_prior_f1) + '\n') 
    # f.write('GENIE normal: ' + str(genie_normal_f1) + '\n')
    # f.write('NI normal: ' + str(ni_normal_f1) + '\n')             
    f.write('\n')  
    f.write('\n\n')
    f.write('AUPR\'\ score:\n')
    f.write('No priors: ' + str(no_prior_score) + '\n')
    f.write('Full priors: ' + str(full_prior_score) + '\n')
    f.write('\n')  
    f.write('\n\n')
    f.write('\nNo priors:\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_no)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_no)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_no)) + '\n')
    f.write('\n')
    f.write('\nFull priors:\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_full)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_full)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_full)) + '\n')
    f.write('\n')  
    # f.write('\nGenie priors:\n')
    # f.write('Agreed upon gene connections = ' + str(len(agreed_g_prior)) + '\n')
    # f.write('\n')
    # f.write('Missed gene connections = ' + str(len(orig_found_g_prior)) + '\n')
    # f.write('\n')
    # f.write('Incorrect gene connections = ' + str(len(sim_found_g_prior)) + '\n')
    # f.write('\n')
    # f.write('\nNlnet priors:\n')
    # f.write('Agreed upon gene connections = ' + str(len(agreed_nlnet)) + '\n')
    # f.write('\n')
    # f.write('Missed gene connections = ' + str(len(orig_found_nlnet)) + '\n')
    # f.write('\n')
    # f.write('Incorrect gene connections = ' + str(len(sim_found_nlnet)) + '\n')
    # f.write('\n')
    # f.write('\nGenie normal\n')
    # f.write('Agreed upon gene connections = ' + str(len(agreed_g_norm)) + '\n')
    # f.write('\n')
    # f.write('Missed gene connections = ' + str(len(orig_found_g_norm)) + '\n')
    # f.write('\n')
    # f.write('Incorrect gene connections = ' + str(len(sim_found_g_norm)) + '\n')
    # f.write('\n')
    # f.write('\nNI normal\n')
    # f.write('Agreed upon gene connections = ' + str(len(agreed_ni)) + '\n')
    # f.write('\n')
    # f.write('Missed gene connections = ' + str(len(orig_found_ni)) + '\n')
    # f.write('\n')
    # f.write('Incorrect gene connections = ' + str(len(sim_found_ni)) + '\n')
    f.write('\n')
    # f.write(str(agreed_no))
    # f.write('\n')
    # f.write(str(len(set(agreed_no))))
    # f.write('\n')
    # f.write(str(agreed_full))
    # f.write('\n')
    # f.write(str(len(agreed_full)))
    # f.write('\n')


