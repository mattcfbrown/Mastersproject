#This file will output metrics needed to write up section 4.3.3 of my report

import numpy as np
import sklearn.metrics as metrics
import sys

#Firstly let's read in the .csv file
no_prior = np.genfromtxt(sys.argv[1], delimiter = ",")
full_prior = np.genfromtxt(sys.argv[2], delimiter = ",")
genie_prior = np.genfromtxt(sys.argv[3], delimiter = ",")
genie = np.genfromtxt(sys.argv[4], delimiter = ",")

#Now the actual matrix
original = np.genfromtxt(sys.argv[5], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))

#Now converying the matricies into a long list where AUROC comparisons can be made
original_list = [x for row in original for x in row]
no_prior_list = [x for row in no_prior for x in row]
full_prior_list = [x for row in full_prior for x in row]
genie_prior_list = [x for row in genie_prior for x in row]
genie_list = [x for row in genie for x in row]

top_idx = np.argsort(genie_list)[-(round(len(genie_list)*0.05)):]
min_val = min([genie_list[i] for i in top_idx])
for j in range(len(genie_list)):
    if genie_list[j] >= min_val:
        genie_list[j] = 1
    else:
        genie_list[j] = 0
for i in range(num_genes):
    for j in range(num_genes):
        if genie[i][j] >= min_val:
            genie[i][j] = 1
        else:
            genie[i][j] = 0

#Getting Matthews correlation:
no_prior_matt = metrics.matthews_corrcoef(original_list,no_prior_list)
full_prior_matt = metrics.matthews_corrcoef(original_list,full_prior_list)
genie_prior_matt = metrics.matthews_corrcoef(original_list,genie_prior_list)
genie_matt = metrics.matthews_corrcoef(original_list,genie_list)

#Get a F1 score
no_prior_f1 = metrics.f1_score(original_list,no_prior_list)
full_prior_f1 = metrics.f1_score(original_list,full_prior_list)
genie_prior_f1 = metrics.f1_score(original_list,genie_prior_list)
genie_f1 = metrics.f1_score(original_list,genie_list)

#Now calculate the Area under the curve
#Note for now this only calculates the NLNET AUROC, but other metrics can be added once this has been completed. 
precision_no_prior, recall_no_prior, thresholds_no_prior = metrics.precision_recall_curve(original_list,no_prior_list)
precision_full_prior, recall_full_prior, thresholds_full_prior = metrics.precision_recall_curve(original_list,full_prior_list)
precision_genie, recall_genie, thresholds_genie = metrics.precision_recall_curve(original_list,genie_prior_list)
precision_nlent, recall_nlnet, thresholds_nlnet = metrics.precision_recall_curve(original_list,genie_list)

no_prior_score = metrics.auc(recall_no_prior,precision_no_prior)
full_prior_score = metrics.auc(recall_full_prior,precision_full_prior)
genie_prior_score = metrics.auc(recall_genie,precision_genie)
genie_score = metrics.auc(recall_nlnet,precision_nlent)

#-----------------------------------Now getting the other stuff--------------------------------

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
agreed_g_norm = [] #Connections both found
sim_found_g_norm = [] #Connections only found by EB
orig_found_g_norm = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if genie_prior [i][j] == 1:
            sim_found_g_norm.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_g_norm.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_g_norm = set(sim_found_g_norm).intersection(set(orig_found_g_norm))
agreed_g_norm = list(agreed_g_norm)
sim_found_g_norm = [x for x in sim_found_g_norm if x not in agreed_g_norm]
orig_found_g_norm = [x for x in orig_found_g_norm if x not in agreed_g_norm]

#Nlnet prior
agreed_nlnet = [] #Connections both found
sim_found_nlnet = [] #Connections only found by EB
orig_found_nlnet = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if genie [i][j] == 1:
            sim_found_nlnet.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_nlnet.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_nlnet = set(sim_found_nlnet).intersection(set(orig_found_nlnet))
agreed_nlnet = list(agreed_nlnet)
sim_found_nlnet = [x for x in sim_found_nlnet if x not in agreed_nlnet]
orig_found_nlnet = [x for x in orig_found_nlnet if x not in agreed_nlnet]

#Now that we have our 6, we can now write it into text files
id = sys.argv[6]
w0 = sys.argv[7]
matthew_name = 'matthew_' + str(id) + '_' + str(w0) + '.txt'
with open(matthew_name, 'w') as f:
    f.write('Matthew\'\ correlation:\n')
    f.write('No priors: ' + str(no_prior_matt) + '\n')
    f.write('Full priors: ' + str(full_prior_matt) + '\n')
    f.write('GENIE3 priors: ' + str(genie_prior_matt) + '\n')
    f.write('GENIE3: ' + str(genie_matt) + '\n')      
    f.write('\n')
    f.write('F1\'\ score:\n')
    f.write('No priors: ' + str(no_prior_f1) + '\n')
    f.write('Full priors: ' + str(full_prior_f1) + '\n')
    f.write('GENIE3 priors: ' + str(genie_prior_f1) + '\n')
    f.write('GENIE3: ' + str(genie_f1) + '\n')          
    f.write('\n')  
    f.write('\n\n')
    f.write('AUPR\'\ score:\n')
    f.write('No priors: ' + str(no_prior_score) + '\n')
    f.write('Full priors: ' + str(full_prior_score) + '\n')
    f.write('GENIE3 priors: ' + str(genie_prior_score) + '\n')
    f.write('GENIE3: ' + str(genie_score) + '\n') 
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
    f.write('\n')
    f.write('\nGENIE3:\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_nlnet)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_nlnet)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_nlnet)) + '\n')
    f.write('\n')
    f.write('\GENIE3 priors\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_g_norm)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_g_norm)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_g_norm)) + '\n')
    f.write('\n')
