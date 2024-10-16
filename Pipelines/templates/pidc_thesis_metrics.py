#Here I will have a simple metrics output thing

#Libraries which I will be using
import numpy as np
import sys
import sklearn.metrics as metrics


#I will read in the five matricies into a list
#I will then for loop it all so I am working on it one at a time
#And then I can write it all into a text file one at a time

#Some test data which I will use to test out my program

#Firstly read in the data
no_priors = np.genfromtxt(sys.argv[1], delimiter = ",")
full_priors = np.genfromtxt(sys.argv[2], delimiter = ",")
genie3_priors = np.genfromtxt(sys.argv[3], delimiter = ",")
genie3 = np.genfromtxt(sys.argv[4], delimiter = ",")
pidc = np.genfromtxt(sys.argv[7], delimiter = ",")
# nlnet_priors = np.genfromtxt(sys.argv[5], delimiter = ",")
#The original
original = np.genfromtxt(sys.argv[5], delimiter = ",")
original[0][0] = 0
num_genes = int(np.size(original, 1))
print(num_genes)

#With genie3, we need to make it so that the 0s and 1s are there
min_val = sorted([x for row in genie3 for x in row],reverse=True)
min_val = min_val[round(num_genes*num_genes*0.05)]
for i in range(num_genes):
    for j in range(num_genes):
        if genie3[i][j] >= min_val or genie3[j][i] == 1:
            genie3[i][j] = 1
            genie3[j][i] = 1
        else:
            genie3[i][j] = 0

#Put everything into a list:
names = ['No priors',
         'Full priors',
         'GENIE3 priors',
         'GENIE3',
         'PIDC']

#List of outputs
outputs = [no_priors,full_priors,genie3_priors,genie3,pidc]

#Things we will be using
mcc = []
f1 = []
AUPR = []

#Now found and missing data
agreed= [] #Connections both found
sim_found = [] #Connections only found by EB
orig_found = [] #Connections not found by EB

#We now create a gene list for this data, to work out the connected edges
gene_list = []
#We create our gene list
for i in range(num_genes):
    gene_list.append('T' + str(i))

#Let's loop it all
for data in outputs:
    data_list = [x for row in data for x in row]
    original_list = [x for row in original for x in row]
    #Matthew's correlation
    mcc.append(metrics.matthews_corrcoef(original_list,data_list))
    #F1 score:
    f1.append(metrics.f1_score(original_list,data_list))
    #AUPR:
    precision, recall, thresholds = metrics.precision_recall_curve(original_list,data_list)
    AUPR.append(metrics.auc(recall,precision))

    #Now the connection information
    agreed_holder = []
    sim_found_holder = []
    orig_found_holder = []
    #Let's get the simulated and original ones
    for i in range(num_genes):
        for j in range(i,num_genes):
            if data [i][j] == 1:
                sim_found_holder.append(gene_list[i] + ' --- ' + gene_list[j])
            if original[i][j] == 1:
                orig_found_holder.append(gene_list[i] + ' --- ' + gene_list[j])     

    #We get agreement
    agreed_holder = set(sim_found_holder).intersection(set(orig_found_holder))
    agreed.append(len(list(agreed_holder)))
    sim_found.append(len([x for x in sim_found_holder if x not in agreed_holder]))
    orig_found.append(len([x for x in orig_found_holder if x not in agreed_holder]))     

#We now read all which is above into a file
file_name = 'metrics_' + sys.argv[6] + '.txt'
with open(file_name, 'w') as f:
    #Scores
    f.write('MCC score: ' + '\n\n')
    for i in range(len(names)):
        f.write(names[i] + ': ')
        f.write(str(mcc[i]) + '\n')
    f.write('\n')
    f.write('F1 score: ' + '\n\n')
    for i in range(len(names)):
        f.write(names[i] + ': ')
        f.write(str(f1[i]) + '\n')
    f.write('\n')
    f.write('AUPR score: ' + '\n\n')
    for i in range(len(names)):
        f.write(names[i] + ': ')
        f.write(str(AUPR[i]) + '\n')
    f.write('\n')
    #Now time to write in the missed
    for i in range(len(names)):
        f.write(names[i] + ' connection information\n\n')
        f.write('Found: ')
        f.write(str(agreed[i]))
        f.write('\n')
        f.write('Incorrect: ')
        f.write(str(sim_found[i]))
        f.write('\n')
        f.write('Missed: ')
        f.write(str(orig_found[i]))
        f.write('\n')
        f.write('\n\n')
        



