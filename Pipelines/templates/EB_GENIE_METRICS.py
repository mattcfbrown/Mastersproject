#This file has three goals:
#1 Get AUPR Scores for the GENIE priors and the Genie data
#2 Get a table with the following data
    #Connections which GENIE and GENIE priors agree on
    #Connections which GENIE found and GENIE priors did not
    #Connections which GENIE did not find but GENIE3 did
#3 Plots, showing the 3 above senarios
    #Show the following data:
        #Data
        #Prior function
        #Posterior due to the data

#--------------------Part 0-------------------------------------
#We load in the libaries we will be using
import sklearn.metrics as metrics
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

#We now load in the data
genie = np.genfromtxt(sys.argv[1], delimiter = ",")                     #This is the output from GENIE
genie_prior = np.genfromtxt(sys.argv[2], delimiter = ",")               #This is the output using GENIE priors
original = np.genfromtxt(sys.argv[3], delimiter = ",")                  #This is the actual network
data = np.genfromtxt(sys.argv[4], delimiter = ",")                      #This is a list of the data
gene_list = np.genfromtxt(sys.argv[5], delimiter = "\t",dtype=None)     #This is a list of all the genes
prior_info = np.genfromtxt(sys.argv[6], delimiter = ",")                #This is the prior information

id  = sys.argv[7]                                                       #This is the id, the way we can tell them apart

#--------------------Part 1-------------------------------------
#We now try to calculate the AUPR scores

#Put the scores into a workable form
genie_list = [x for row in genie for x in row]
genie_prior_list = [x for row in genie_prior for x in row]
original_list = [x for row in original for x in row]

#We calculate the AUPR scores
original_list[0] = 0.0
precision_genie, recall_genie, thresholds_genie = metrics.precision_recall_curve(original_list,genie_list)
precision_prior, recall_prior, thresholds_prior = metrics.precision_recall_curve(original_list,genie_prior_list)


#--------------------Part 2----------------------------------
#Here we now create a list of genes which have agreement and disagreement

#We have three lists
agreed = []                 #Correct edges found by both GENIE3 and GENIE3 priors
genie_found = []            #Connections only found by GENIE3
prior_found = []            #Connections not found by GENIE3 priors
incorrect_genie = []        #Incorrectly found genie
incorrect_prior = []        #Incorrectly found prior

for i in range(genie.shape[0]):
    for j in range(i,genie.shape[1]):
        #Agreement between both
        if genie[i][j] == genie_prior[i][j] and original[i][j] == 1:
            agreed.append('T' + str(i) + '----' + 'T' + str(j))
        #GENIE3 only
        elif genie[i][j] == 1 and genie_prior[i][j] == 0 and original[i][j] == 1:
            genie_found.append('T' + str(i) + '----' + 'T' + str(j))
        #GENIE3 priors only
        elif genie[i][j] == 0 and genie_prior[i][j] == 1 and original[i][j] == 1:
            prior_found.append('T' + str(i) + '----' + 'T' + str(j))
        if genie[i][j] == 1 and original[i][j] == 0:
            incorrect_genie.append('T' + str(i) + '----' + 'T' + str(j))
        if genie_prior[i][j] == 1 and original[i][j] == 0:
            incorrect_prior.append('T' + str(i) + '----' + 'T' + str(j))            


#We now write this into a file
file_name = 'metrics_' + str(id) +  '.txt'
with open(file_name, 'w') as f:
    f.write('AUPR SCORES:\n\n')
    f.write('AUPR score (GENIE): ' + str(metrics.auc(recall_genie,precision_genie)) + '\n')
    f.write('AUPR score (Priors): ' + str(metrics.auc(recall_prior,precision_prior)) + '\n')
    f.write('\n')
    f.write('--------------------------------------------------------------------------\n')   
    f.write('Agreed upon gene connections (' + str(len(agreed)) + '):  \n')
    for connection in agreed:
        f.write(connection + '\n')
    f.write('\n')
    f.write('GENIE connections (' + str(len(genie_found)) + '):  \n')
    for connection in genie_found:
        f.write(connection + '\n')
    f.write('\n')
    f.write('PRIOR connections (' + str(len(prior_found)) + '):  \n')
    for connection in prior_found:
        f.write(connection + '\n')
    f.write('\n')
    f.write('--------------------------------------------------------------------------\n')
    f.write('Incorrect Genie connections (' + str(len(incorrect_genie)) + '):  \n')
    for connection in incorrect_genie:
        f.write(connection + '\n')
    f.write('\n')
    f.write('Incorrect Prior connections (' + str(len(incorrect_prior)) + '):  \n')
    for connection in incorrect_prior:
        f.write(connection + '\n')
    f.write('\n')


#--------------------Part 3----------------------------------
#This is where the plots are occuring
#We will have a 2x2 plot, which display the following:
    #Top left     = Prior function
    #Bottom left  = Prior value distribution
    #Top right    = Data (The test statistics found)
    #Bottom right = Posterior distribution

#Firstly classify each point in the gene list as 1-5, based on the orderings above
#we also get their prior value
a_priors = []
g_priors = []
p_priors = []
i_g_priors = []
i_p_priors = []
connection_type = []
for genes in gene_list:
    gene1 = genes[0].decode('UTF-8')
    gene2 = genes[1].decode('UTF-8')
    search = gene1 + '----' + gene2
    if search in agreed:
        connection_type.append(1)
        a_priors.append(max(prior_info[int(gene1[1:]),int(gene2[1:])],
                            prior_info[int(gene2[1:]),int(gene1[1:])]) *2.2)
    elif search in genie_found:
        connection_type.append(2)
        g_priors.append(max(prior_info[int(gene1[1:]),int(gene2[1:])],
                            prior_info[int(gene2[1:]),int(gene1[1:])]) *2.2)
    elif search in prior_found:
        connection_type.append(3)
        p_priors.append(max(prior_info[int(gene1[1:]),int(gene2[1:])],
                            prior_info[int(gene2[1:]),int(gene1[1:])]) *2.2)
    elif search in incorrect_genie:
        connection_type.append(4)
        i_g_priors.append(max(prior_info[int(gene1[1:]),int(gene2[1:])],
                            prior_info[int(gene2[1:]),int(gene1[1:])]) *2.2)
    elif search in incorrect_prior:
        connection_type.append(5)
        i_p_priors.append(max(prior_info[int(gene1[1:]),int(gene2[1:])],
                            prior_info[int(gene2[1:]),int(gene1[1:])]) *2.2)
    else:
        connection_type.append(0)

#Write out the prior function:
def prior_fn(x,w0,multi):
    if x == 0.0:
        return math.exp(w0) / ( math.exp(w0) + math.exp(x) )
    elif 0.0*2.2 < x <0.1*2.2:
        return math.exp(w0) / ( math.exp(w0) + 0*math.exp(x) )
    elif 0.1*2.2 <= x < 0.5*2.2:
        return math.exp(w0) / ( math.exp(w0) + multi*0.5*math.exp(x) )        
    else:
       return  math.exp(w0) / ( math.exp(w0) + multi*math.exp(x) )  

#Get their index locations
agreed_index = [i for i, value in enumerate(connection_type) if value == 1]
genie_index = [i for i, value in enumerate(connection_type) if value == 2]
prior_index = [i for i, value in enumerate(connection_type) if value == 3]
incorrect_genie_index = [i for i, value in enumerate(connection_type) if value == 4]
incorrect_prior_index = [i for i, value in enumerate(connection_type) if value == 5]

#We will firstly test it works for the agreed as a proof of concept, then expand for others
#This gets the prior functions
x = np.linspace(0, 2.2, 300)
function_val = []
for x_val in x:
    function_val.append(prior_fn(x_val,2.2,10.0))

#----------------------------------Part 4: Agreed-------------------------------------------
fig, ax = plt.subplots(2, 2)
ax[0,0].plot(x,function_val)                                                  #Prior function
ax[0,0].set_xlabel('x')
ax[0,0].set_ylabel('Alpha value')
ax[0,0].set_title('Prior function')

colours_agreed = cm.rainbow(np.linspace(0, 1, len(agreed_index)))
y_agreed = np.random.rand(len(agreed_index))/2

ax[0,1].scatter(data[2][agreed_index], y_agreed,                                             #Data Linegraph
                c=colours_agreed,
                s = 10)
ax[0,1].set_xlabel('Test statistic values')
ax[0,1].set_ylabel(' ')
ax[0,1].set_title('Line graph of test statistics')
ax[0,1].set_ylim([0, 0.75])


for i, txt in enumerate([value for value in range(1,len(agreed_index)+1,1)]):
    ax[0,1].annotate(txt, (data[2][agreed_index][i], y_agreed[i]),fontsize=5)

ax[1,0].scatter(a_priors,y_agreed,                                                               #Prior Linegraph
                c=colours_agreed,
                s = 10)
ax[1,0].set_xlabel('Prior values')
ax[1,0].set_ylabel(' ')
ax[1,0].set_title('Line graph of Prior values')

ax[1,0].set_xlim([0, 2.2])
ax[1,0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(agreed_index)+1,1)]):
    ax[1,0].annotate(txt, (a_priors[i], y_agreed[i]),fontsize=5)

ax[1,1].scatter(1 - (a_priors)*data[0][agreed_index]/data[1][agreed_index], y_agreed,            #Posterior Linegraph
                c=colours_agreed,
                s = 10)
ax[1,1].set_xlabel('Posterior')
ax[1,1].set_ylabel(' ')
ax[1,1].set_title('Line graph of the Posterior')

ax[1,1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(agreed_index)+1,1)]):
    ax[1,1].annotate(txt, ((1 - (a_priors)*data[0][agreed_index]/data[1][agreed_index])[i], y_agreed[i]),
                     fontsize=5)

fig.suptitle('Plots for agreed upon connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Agreed_" + str(id) + ".pdf", format="pdf")


#----------------------------------Part 4: GENIE-------------------------------------------
fig, ax = plt.subplots(2, 2)
ax[0,0].plot(x,function_val)                                                  #Prior function
ax[0,0].set_xlabel('x')
ax[0,0].set_ylabel('Alpha value')
ax[0,0].set_title('Prior function')

colours_genie = cm.rainbow(np.linspace(0, 1, len(genie_index)))
y_genie = np.random.rand(len(genie_index))/2

ax[0,1].scatter(data[2][genie_index], y_genie,                                             #Data Linegraph
                c=colours_genie,
                s = 10)
ax[0,1].set_xlabel('Test statistic values')
ax[0,1].set_ylabel(' ')
ax[0,1].set_title('Line graph of test statistics')
ax[0,1].set_ylim([0, 0.75])


for i, txt in enumerate([value for value in range(1,len(genie_index)+1,1)]):
    ax[0,1].annotate(txt, (data[2][genie_index][i], y_genie[i]),fontsize=5)

ax[1,0].scatter(g_priors,y_genie,                                                               #Prior Linegraph
                c=colours_genie,
                s = 10)
ax[1,0].set_xlabel('Prior values')
ax[1,0].set_ylabel(' ')
ax[1,0].set_title('Line graph of Prior values')

ax[1,0].set_xlim([0, 2.2])
ax[1,0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(genie_index)+1,1)]):
    ax[1,0].annotate(txt, (g_priors[i], y_genie[i]),fontsize=5)

ax[1,1].scatter(1 - (g_priors)*data[0][genie_index]/data[1][genie_index], y_genie,            #Posterior Linegraph
                c=colours_genie,
                s = 10)
ax[1,1].set_xlabel('Posterior')
ax[1,1].set_ylabel(' ')
ax[1,1].set_title('Line graph of the Posterior')

ax[1,1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(genie_index)+1,1)]):
    ax[1,1].annotate(txt, ((1 - (g_priors)*data[0][genie_index]/data[1][genie_index])[i], y_genie[i]),
                     fontsize=5)

fig.suptitle('Plots for Genie only connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Genie_" + str(id) + ".pdf", format="pdf")


#----------------------------------Part 4: Prior-------------------------------------------
fig, ax = plt.subplots(2, 2)
ax[0,0].plot(x,function_val)                                                  #Prior function
ax[0,0].set_xlabel('x')
ax[0,0].set_ylabel('Alpha value')
ax[0,0].set_title('Prior function')

colours_prior = cm.rainbow(np.linspace(0, 1, len(prior_index)))
y_prior = np.random.rand(len(prior_index))/2

ax[0,1].scatter(data[2][prior_index], y_prior,                                             #Data Linegraph
                c=colours_prior,
                s = 10)
ax[0,1].set_xlabel('Test statistic values')
ax[0,1].set_ylabel(' ')
ax[0,1].set_title('Line graph of test statistics')
ax[0,1].set_ylim([0, 0.75])


for i, txt in enumerate([value for value in range(1,len(prior_index)+1,1)]):
    ax[0,1].annotate(txt, (data[2][prior_index][i], y_prior[i]),fontsize=5)

ax[1,0].scatter(p_priors,y_prior,                                                               #Prior Linegraph
                c=colours_prior,
                s = 10)
ax[1,0].set_xlabel('Prior values')
ax[1,0].set_ylabel(' ')
ax[1,0].set_title('Line graph of Prior values')

ax[1,0].set_xlim([0, 2.2])
ax[1,0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(prior_index)+1,1)]):
    ax[1,0].annotate(txt, (p_priors[i], y_prior[i]),fontsize=5)

ax[1,1].scatter(1 - (p_priors)*data[0][prior_index]/data[1][prior_index], y_prior,            #Posterior Linegraph
                c=colours_prior,
                s = 10)
ax[1,1].set_xlabel('Posterior')
ax[1,1].set_ylabel(' ')
ax[1,1].set_title('Line graph of the Posterior')

ax[1,1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(prior_index)+1,1)]):
    ax[1,1].annotate(txt, ((1 - (p_priors)*data[0][prior_index]/data[1][prior_index])[i], y_prior[i]),
                     fontsize=5)

fig.suptitle('Plots for Prior only connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Prior_" + str(id) + ".pdf", format="pdf")


#----------------------------------Part 4: Incorrect Genie----------------------------------
fig, ax = plt.subplots(2, 2)
ax[0,0].plot(x,function_val)                                                  #Prior function
ax[0,0].set_xlabel('x')
ax[0,0].set_ylabel('Alpha value')
ax[0,0].set_title('Prior function')

colours_incorrect_genie = cm.rainbow(np.linspace(0, 1, len(incorrect_genie_index)))
y_incorrect_genie = np.random.rand(len(incorrect_genie_index))/2

ax[0,1].scatter(data[2][incorrect_genie_index], y_incorrect_genie,                                             #Data Linegraph
                c=colours_incorrect_genie,
                s = 10)
ax[0,1].set_xlabel('Test statistic values')
ax[0,1].set_ylabel(' ')
ax[0,1].set_title('Line graph of test statistics')
ax[0,1].set_ylim([0, 0.75])


for i, txt in enumerate([value for value in range(1,len(incorrect_genie_index)+1,1)]):
    ax[0,1].annotate(txt, (data[2][incorrect_genie_index][i], y_incorrect_genie[i]),fontsize=5)

ax[1,0].scatter(i_g_priors,y_incorrect_genie,                                                               #Prior Linegraph
                c=colours_incorrect_genie,
                s = 10)
ax[1,0].set_xlabel('Prior values')
ax[1,0].set_ylabel(' ')
ax[1,0].set_title('Line graph of Prior values')

ax[1,0].set_xlim([0, 2.2])
ax[1,0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_genie_index)+1,1)]):
    ax[1,0].annotate(txt, (i_g_priors[i], y_incorrect_genie[i]),fontsize=5)

ax[1,1].scatter(1 - (i_g_priors)*data[0][incorrect_genie_index]/data[1][incorrect_genie_index], y_incorrect_genie,            #Posterior Linegraph
                c=colours_incorrect_genie,
                s = 10)
ax[1,1].set_xlabel('Posterior')
ax[1,1].set_ylabel(' ')
ax[1,1].set_title('Line graph of the Posterior')

ax[1,1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_genie_index)+1,1)]):
    ax[1,1].annotate(txt, ((1 - (i_g_priors)*data[0][incorrect_genie_index]/data[1][incorrect_genie_index])[i], y_incorrect_genie[i]),
                     fontsize=5)

fig.suptitle('Plots for Incorrect Genie connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Incorrect_Genie_" + str(id) + ".pdf", format="pdf")


#----------------------------------Part 4: Incorrect Prior----------------------------------
fig, ax = plt.subplots(2, 2)
ax[0,0].plot(x,function_val)                                                            #Prior function
ax[0,0].set_xlabel('x')
ax[0,0].set_ylabel('Alpha value')
ax[0,0].set_title('Prior function')

colours_incorrect_prior = cm.rainbow(np.linspace(0, 1, len(incorrect_prior_index)))
y_incorrect_prior = np.random.rand(len(incorrect_prior_index))/2

ax[0,1].scatter(data[2][incorrect_prior_index], y_incorrect_prior,                                             #Data Linegraph
                c=colours_incorrect_prior,
                s = 10)
ax[0,1].set_xlabel('Test statistic values')
ax[0,1].set_ylabel(' ')
ax[0,1].set_title('Line graph of test statistics')
ax[0,1].set_ylim([0, 0.75])


for i, txt in enumerate([value for value in range(1,len(incorrect_prior_index)+1,1)]):
    ax[0,1].annotate(txt, (data[2][incorrect_prior_index][i], y_incorrect_prior[i]),fontsize=5)

ax[1,0].scatter(i_p_priors,y_incorrect_prior,                                                               #Prior Linegraph
                c=colours_incorrect_prior,
                s = 10)
ax[1,0].set_xlabel('Prior values')
ax[1,0].set_ylabel(' ')
ax[1,0].set_title('Line graph of Prior values')

ax[1,0].set_xlim([0, 2.2])
ax[1,0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_prior_index)+1,1)]):
    ax[1,0].annotate(txt, (i_p_priors[i], y_incorrect_prior[i]),fontsize=5)

ax[1,1].scatter(1 - (i_p_priors)*data[0][incorrect_prior_index]/data[1][incorrect_prior_index], y_incorrect_prior,            #Posterior Linegraph
                c=colours_incorrect_prior,
                s = 10)
ax[1,1].set_xlabel('Posterior')
ax[1,1].set_ylabel(' ')
ax[1,1].set_title('Line graph of the Posterior')

ax[1,1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_prior_index)+1,1)]):
    ax[1,1].annotate(txt, ((1 - (i_p_priors)*data[0][incorrect_prior_index]/data[1][incorrect_prior_index])[i], y_incorrect_prior[i]),
                     fontsize=5)

fig.suptitle('Plots for Incorrect Genie connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Incorrect_Prior_" + str(id) + ".pdf", format="pdf")

#-------------------Null and Mixture densities----------------------
#Agreed:


fig, ax = plt.subplots(1, 3)
ax[0].scatter(data[0][agreed_index],y_agreed,
            c=colours_agreed,
            s = 10)
ax[0].set_xlabel('Null values')
ax[0].set_ylabel(' ')
ax[0].set_title('Null distribution')
ax[0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(agreed_index)+1,1)]):
    ax[0].annotate(txt, (data[0][agreed_index][i], y_agreed[i]),
                     fontsize=5)


ax[1].scatter(data[1][agreed_index],y_agreed,
            c=colours_agreed,
            s = 10)
ax[1].set_xlabel('Mixture values')
ax[1].set_ylabel(' ')
ax[1].set_title('Mixture distribution')
ax[1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(agreed_index)+1,1)]):
    ax[1].annotate(txt, (data[1][agreed_index][i], y_agreed[i]),
                     fontsize=5)

ax[2].scatter(data[0][agreed_index]/data[1][agreed_index],y_agreed,
            c=colours_agreed,
            s = 10)
ax[2].set_xlabel('Null/Mixture values')
ax[2].set_ylabel(' ')
ax[2].set_title('Null/Mixture distribution')
ax[2].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(agreed_index)+1,1)]):
    ax[2].annotate(txt, ((data[0][agreed_index]/data[1][agreed_index])[i], y_agreed[i]),
                     fontsize=5)

fig.suptitle('Plots for Agreed connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Null_Mixture_Agreed_" + str(id) + ".pdf", format="pdf")

#Genie:
fig, ax = plt.subplots(1, 3)
ax[0].scatter(data[0][genie_index],y_genie,
            c=colours_genie,
            s = 10)
ax[0].set_xlabel('Null values')
ax[0].set_ylabel(' ')
ax[0].set_title('Null distribution')
ax[0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(genie_index)+1,1)]):
    ax[0].annotate(txt, (data[0][genie_index][i], y_genie[i]),
                     fontsize=5)


ax[1].scatter(data[1][genie_index],y_genie,
            c=colours_genie,
            s = 10)
ax[1].set_xlabel('Mixture values')
ax[1].set_ylabel(' ')
ax[1].set_title('Mixture distribution')
ax[1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(genie_index)+1,1)]):
    ax[1].annotate(txt, (data[1][genie_index][i], y_genie[i]),
                     fontsize=5)

ax[2].scatter(data[0][genie_index]/data[1][genie_index],y_genie,
            c=colours_genie,
            s = 10)
ax[2].set_xlabel('Null/Mixture values')
ax[2].set_ylabel(' ')
ax[2].set_title('Null/Mixture distribution')
ax[2].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(genie_index)+1,1)]):
    ax[2].annotate(txt, ((data[0][genie_index]/data[1][genie_index])[i], y_genie[i]),
                     fontsize=5)

fig.suptitle('Plots for Genie connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Null_Mixture_Genie_" + str(id) + ".pdf", format="pdf")

#prior:
fig, ax = plt.subplots(1, 3)
ax[0].scatter(data[0][prior_index],y_prior,
            c=colours_prior,
            s = 10)
ax[0].set_xlabel('Null values')
ax[0].set_ylabel(' ')
ax[0].set_title('Null distribution')
ax[0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(prior_index)+1,1)]):
    ax[0].annotate(txt, (data[0][prior_index][i], y_prior[i]),
                     fontsize=5)


ax[1].scatter(data[1][prior_index],y_prior,
            c=colours_prior,
            s = 10)
ax[1].set_xlabel('Mixture values')
ax[1].set_ylabel(' ')
ax[1].set_title('Mixture distribution')
ax[1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(prior_index)+1,1)]):
    ax[1].annotate(txt, (data[1][prior_index][i], y_prior[i]),
                     fontsize=5)

ax[2].scatter(data[0][prior_index]/data[1][prior_index],y_prior,
            c=colours_prior,
            s = 10)
ax[2].set_xlabel('Null/Mixture values')
ax[2].set_ylabel(' ')
ax[2].set_title('Null/Mixture distribution')
ax[2].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(prior_index)+1,1)]):
    ax[2].annotate(txt, ((data[0][prior_index]/data[1][prior_index])[i], y_prior[i]),
                     fontsize=5)

fig.suptitle('Plots for prior connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Null_Mixture_Prior_" + str(id) + ".pdf", format="pdf")

#Incorrect Genie:
fig, ax = plt.subplots(1, 3)
ax[0].scatter(data[0][incorrect_genie_index],y_incorrect_genie,
            c=colours_incorrect_genie,
            s = 10)
ax[0].set_xlabel('Null values')
ax[0].set_ylabel(' ')
ax[0].set_title('Null distribution')
ax[0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_genie_index)+1,1)]):
    ax[0].annotate(txt, (data[0][incorrect_genie_index][i], y_incorrect_genie[i]),
                     fontsize=5)


ax[1].scatter(data[1][incorrect_genie_index],y_incorrect_genie,
            c=colours_incorrect_genie,
            s = 10)
ax[1].set_xlabel('Mixture values')
ax[1].set_ylabel(' ')
ax[1].set_title('Mixture distribution')
ax[1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_genie_index)+1,1)]):
    ax[1].annotate(txt, (data[1][incorrect_genie_index][i], y_incorrect_genie[i]),
                     fontsize=5)

ax[2].scatter(data[0][incorrect_genie_index]/data[1][incorrect_genie_index],y_incorrect_genie,
            c=colours_incorrect_genie,
            s = 10)
ax[2].set_xlabel('Null/Mixture values')
ax[2].set_ylabel(' ')
ax[2].set_title('Null/Mixture distribution')
ax[2].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_genie_index)+1,1)]):
    ax[2].annotate(txt, ((data[0][incorrect_genie_index]/data[1][incorrect_genie_index])[i], y_incorrect_genie[i]),
                     fontsize=5)

fig.suptitle('Plots for Incorrect Genie connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Null_Mixture_Incorrect_Genie_" + str(id) + ".pdf", format="pdf")

#Incorrect Genie:
fig, ax = plt.subplots(1, 3)
ax[0].scatter(data[0][incorrect_prior_index],y_incorrect_prior,
            c=colours_incorrect_prior,
            s = 10)
ax[0].set_xlabel('Null values')
ax[0].set_ylabel(' ')
ax[0].set_title('Null distribution')
ax[0].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_prior_index)+1,1)]):
    ax[0].annotate(txt, (data[0][incorrect_prior_index][i], y_incorrect_prior[i]),
                     fontsize=5)


ax[1].scatter(data[1][incorrect_prior_index],y_incorrect_prior,
            c=colours_incorrect_prior,
            s = 10)
ax[1].set_xlabel('Mixture values')
ax[1].set_ylabel(' ')
ax[1].set_title('Mixture distribution')
ax[1].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_prior_index)+1,1)]):
    ax[1].annotate(txt, (data[1][incorrect_prior_index][i], y_incorrect_prior[i]),
                     fontsize=5)

ax[2].scatter(data[0][incorrect_prior_index]/data[1][incorrect_prior_index],y_incorrect_prior,
            c=colours_incorrect_prior,
            s = 10)
ax[2].set_xlabel('Null/Mixture values')
ax[2].set_ylabel(' ')
ax[2].set_title('Null/Mixture distribution')
ax[2].set_ylim([0, 0.75])

for i, txt in enumerate([value for value in range(1,len(incorrect_prior_index)+1,1)]):
    ax[2].annotate(txt, ((data[0][incorrect_prior_index]/data[1][incorrect_prior_index])[i], y_incorrect_prior[i]),
                     fontsize=5)


fig.suptitle('Plots for Incorrect Prior connections for ' + str(id) + ' cells')
plt.tight_layout()

#We save the plot here
plt.savefig("Null_Mixture_Incorrect_prior_" + str(id) + ".pdf", format="pdf")