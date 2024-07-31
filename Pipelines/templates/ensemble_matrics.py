#This file here gets the metrics from  the ensemble method I have developed.

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import sys

#-----------------------------Read in the data-----------------------------------
original = np.genfromtxt(sys.argv[1], delimiter = ",")
original[0][0] = 0
eb = np.genfromtxt(sys.argv[2], delimiter = ",")
genie = np.genfromtxt(sys.argv[3], delimiter = ",")
mi = np.genfromtxt(sys.argv[4], delimiter = ",")
puc = np.genfromtxt(sys.argv[5], delimiter = ",")
clr = np.genfromtxt(sys.argv[6], delimiter = ",")
ensemble = np.genfromtxt(sys.argv[7], delimiter = ",")

#-----------------------------Convert it into a form we want to us---------------
original_list = [x for row in original for x in row]
eb_list = [x for row in eb for x in row]
genie_list = [x for row in genie for x in row]
mi_list = [x for row in mi for x in row]
puc_list = [x for row in puc for x in row]
clr_list = [x for row in clr for x in row]
ensemble_list = [x for row in ensemble for x in row]

#--------------------------AUPR scores-------------------------------------------
precision_eb, recall_eb, thresholds_eb = metrics.precision_recall_curve(original_list,eb_list)
precision_genie, recall_genie, thresholds_genie = metrics.precision_recall_curve(original_list,genie_list)
precision_mi, recall_mi, thresholds_mi = metrics.precision_recall_curve(original_list,mi_list)
precision_puc, recall_puc, thresholds_puc = metrics.precision_recall_curve(original_list,puc_list)
precision_clr, recall_clr, thresholds_clr = metrics.precision_recall_curve(original_list,clr_list)
precision_ens, recall_ens, thresholds_ens = metrics.precision_recall_curve(original_list,ensemble_list)


eb_score = metrics.auc(recall_eb,precision_eb)
genie_score = metrics.auc(recall_genie,precision_genie)
mi_score = metrics.auc(recall_mi,precision_mi)
puc_score = metrics.auc(recall_puc,precision_puc)
clr_score = metrics.auc(recall_clr,precision_clr)
ens_score = metrics.auc(recall_ens,precision_ens)

#--------------Finding agreed, disagreed, and missed connections------------------
#We now create a gene list for this data, to work out the connected edges
gene_list = []
num_genes = 25
#We create our gene list
for i in range(num_genes):
    gene_list.append('T' + str(i))

#We now write code for each situation
#Empirical Bayes
agreed_eb = [] #Connections both found
sim_found_eb = [] #Connections only found by EB
orig_found_eb = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if eb [i][j] == 1:
            sim_found_eb.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_eb.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_eb = set(sim_found_eb).intersection(set(orig_found_eb))
agreed_eb = list(agreed_eb)
sim_found_eb = [x for x in sim_found_eb if x not in agreed_eb]
orig_found_eb = [x for x in orig_found_eb if x not in agreed_eb]

#GENIE (GENIE = FULL)
agreed_full = [] #Connections both found
sim_found_full = [] #Connections only found by EB
orig_found_full = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if genie[i][j] == 1:
            sim_found_full.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_full.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_full = set(sim_found_full).intersection(set(orig_found_full))
agreed_full = list(agreed_full)
sim_found_full = [x for x in sim_found_full if x not in agreed_full]
orig_found_full = [x for x in orig_found_full if x not in agreed_full]

#MI (MI = GENIE)
agreed_g_norm = [] #Connections both found
sim_found_g_norm = [] #Connections only found by EB
orig_found_g_norm = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if mi[i][j] == 1:
            sim_found_g_norm.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_g_norm.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_g_norm = set(sim_found_g_norm).intersection(set(orig_found_g_norm))
agreed_g_norm = list(agreed_g_norm)
sim_found_g_norm = [x for x in sim_found_g_norm if x not in agreed_g_norm]
orig_found_g_norm = [x for x in orig_found_g_norm if x not in agreed_g_norm]

#PUC (PUC = NLNET)
agreed_nlnet = [] #Connections both found
sim_found_nlnet = [] #Connections only found by EB
orig_found_nlnet = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if puc[i][j] == 1:
            sim_found_nlnet.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_nlnet.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_nlnet = set(sim_found_nlnet).intersection(set(orig_found_nlnet))
agreed_nlnet = list(agreed_nlnet)
sim_found_nlnet = [x for x in sim_found_nlnet if x not in agreed_nlnet]
orig_found_nlnet = [x for x in orig_found_nlnet if x not in agreed_nlnet]

#CLR (CLR = GENIE PRIOR)
agreed_g_prior = [] #Connections both found
sim_found_g_prior = [] #Connections only found by EB
orig_found_g_prior = [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if clr[i][j] == 1:
            sim_found_g_prior.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_g_prior.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_g_prior = set(sim_found_g_prior).intersection(set(orig_found_g_prior))
agreed_g_prior = list(agreed_g_prior)
sim_found_g_prior = [x for x in sim_found_g_prior if x not in agreed_g_prior]
orig_found_g_prior = [x for x in orig_found_g_prior if x not in agreed_g_prior]

#Ensemble
agreed_ens = [] #Connections both found
sim_found_ens = [] #Connections only found by EB
orig_found_ens= [] #Connections not found by EB

#Let's get the simulated and original ones
for i in range(num_genes):
    for j in range(i,num_genes):
        if ensemble[i][j] == 1:
            sim_found_ens.append(gene_list[i] + ' --- ' + gene_list[j])
        if original[i][j] == 1:
            orig_found_ens.append(gene_list[i] + ' --- ' + gene_list[j])            

#We get agreement
agreed_ens = set(sim_found_ens).intersection(set(orig_found_ens))
agreed_ens = list(agreed_ens)
sim_found_ens = [x for x in sim_found_ens if x not in agreed_ens]
orig_found_ens = [x for x in orig_found_ens if x not in agreed_ens]

#------------------------Put it in a text document------------------------------
id = sys.argv[8]
matthew_name = 'matthew_' + str(id) + '.txt'

with open(matthew_name, 'w') as f:
    f.write('AUPR\'\ scores:\n')
    f.write('Empirical Bayes: ' + str(eb_score) + '\n')
    f.write('GENIE: ' + str(genie_score) + '\n')
    f.write('TIGRESS: ' + str(mi_score) + '\n')
    f.write('Pearson: ' + str(puc_score) + '\n') 
    f.write('CLR: ' + str(clr_score) + '\n')   
    f.write('Ensemble: ' + str(ens_score) + '\n')         
    f.write('\n')  
    f.write('\n\n')
    f.write('\nEmprical Bayes:\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_eb)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_eb)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_eb)) + '\n')
    f.write('\n')
    f.write('\nGenie:\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_full)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_full)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_full)) + '\n')
    f.write('\n')  
    f.write('\nTIGRESS:\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_g_prior)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_g_prior)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_g_prior)) + '\n')
    f.write('\n')
    f.write('\Pearson:\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_nlnet)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_nlnet)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_nlnet)) + '\n')
    f.write('\n')
    f.write('\nCLR\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_g_norm)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_g_norm)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_g_norm)) + '\n')
    f.write('\n')
    f.write('Ensemble\n')
    f.write('Agreed upon gene connections = ' + str(len(agreed_ens)) + '\n')
    f.write('\n')
    f.write('Missed gene connections = ' + str(len(orig_found_ens)) + '\n')
    f.write('\n')
    f.write('Incorrect gene connections = ' + str(len(sim_found_ens)) + '\n')
    f.write('\n')