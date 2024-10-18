#This here makes a figure for section 4.5
#I am to create a heatmap for precision and recall (so two in total)

#We firstly need the following libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import seaborn as sns
from statistics import mean
import pandas as pd
import math
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join

path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/metrics'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#We get the names of the methods we wish to test:
method_names = []
file = open(path + '/' + onlyfiles[1], 'r')
for line in file.readlines():
    if 'connection' in line:
        break
    if len(line.rstrip()) > 0 and 'GRNVBEM' not in line and 'AUPR' not in line:
       method_names.append(line.split(':', 1)[0])


#We now want to get the data from each file
network_names = []                       #Holds the names of the networks
AUPR = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the AUPR scores
found_array = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the found scores
missed_array = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the missed scores
incorrect_array = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the incorrect scores
i = 0
for method in onlyfiles:
    #We firstly get all the names
    name = method.replace("Metrics_", "")
    name = name.replace(".txt", "")
    network_names.append(name)
    
    #We are now getting all the values needed
    j = 0
    AUPR_done = 0
    found = 0
    missed = 0
    incorrect = 0
    file = open(path + '/' + method, 'r')
    for line in file.readlines():
        if 'connection' in line:
            AUPR_done = 1
        if AUPR_done == 0:
            if len(line.rstrip()) > 0 and 'GRNVBEM' not in line and 'AUPR' not in line:
                value = line.replace(method_names[j] + ': ', '')
                value = value.replace('\n', '')
                AUPR[i][j] = value
                j = j+1
        else:
            #Get the method
            if 'connection' in line and 'GRNVBEM' not in line:
                grnvbem = 0
                current_method = line.replace(' connection information\n', '')
                j = method_names.index(current_method)
            if 'GRNVBEM' in line:
                incorrect = 0
                found = 0
                missed = 0
                grnvbem = 1
            #Defines where we are
            if 'FOUND' in line and grnvbem == 0:
                found = 1
                incorrect = 0
            elif 'MISSED' in line and grnvbem == 0:
                missed = 1
                found = 0
            elif 'INCORRECT' in line and grnvbem == 0:
                incorrect = 1
                missed = 0
            #We now do the counting
            if found == 1 and '---' in line:
                found_array[i][j] = found_array[i][j] + 1
            elif missed == 1 and '---' in line:
                missed_array[i][j] = missed_array[i][j] + 1
            elif incorrect == 1 and '---' in line:
                incorrect_array[i][j] = incorrect_array[i][j] + 1
    
    i = i+1

#----------------------------Heat map time----------------------

ranking_high = len(method_names)
heat_map_precision = np.zeros((len(onlyfiles),
                     ranking_high)) #Will hold the heat_map values 
general_precision = np.zeros((len(onlyfiles),
                     ranking_high))
normalised_precision = np.zeros((len(onlyfiles),
                     ranking_high))

j = 0
for i in range(len(found_array)):
    TP = found_array[i]
    FP = incorrect_array[i]
    precision = [a/(a+b) for a,b in zip(TP,FP)]
    min_max_precision = precision

    #This is done if using the rank system
    precision = [-1 if math.isnan(x) else x for x in precision]
    array_holder = np.array(precision)
    rank = list(np.argsort((-1*array_holder).argsort()))

    #Editing the rank value
    loc_before = rank.index(0)
    # print(precision[loc_before])
    for i in range(1,ranking_high):
        loc = rank.index(i)
        #If we have the same value, update the rankings
        if precision[loc] == precision[loc_before]:
            rank[loc] = rank[loc_before]
        if precision[loc] == -1:
            rank[loc] = 100
        loc_before = loc

    #This will be done for min_max
    min_precision = min([x for x in min_max_precision if str(x) != 'nan'])
    max_precision = max([x for x in min_max_precision if str(x) != 'nan'])
    for i in range(len(min_max_precision)):
        if str(min_max_precision[i]) == 'nan':
            min_max_precision[i] = 100
        else:
            min_max_precision[i] = (min_max_precision[i]-min_precision)/(max_precision-min_precision)
    
    normalised_precision[j][:] = min_max_precision
    general_precision[j][:] = precision
    heat_map_precision[j][:] = rank
    j = j+1

heat_map_recall = np.zeros((len(onlyfiles),
                     ranking_high)) #Will hold the heat_map values 
general_recall = np.zeros((len(onlyfiles),
                     ranking_high))
normalised_recall = np.zeros((len(onlyfiles),
                     ranking_high))



j = 0
for i in range(len(found_array)):
    TP = found_array[i]
    FN = missed_array[i]
    recall = [a/(a+b) for a,b in zip(TP,FN)]
    for k in range(len(incorrect_array[i])):
        if (incorrect_array[i][k] + found_array[i][k] == 0):
            recall[k] = 100
    min_max_recall = recall
    array_holder = np.array(recall)
    rank = list(np.argsort((-1*array_holder).argsort()))

    #Editing the rank value
    loc_before = rank.index(0)
    for i in range(1,ranking_high):
        loc = rank.index(i)
        #If we have the same value, update the rankings
        if recall[loc] == recall[loc_before]:
            rank[loc] = rank[loc_before]
        loc_before = loc

    print(recall)
    #This will be done for min_max
    min_recall = min([x for x in min_max_recall if x != 100])
    max_recall = max([x for x in min_max_recall if x != 100])
    for i in range(len(min_max_recall)):
        if min_max_recall[i] == 100:
            min_max_recall[i] = 100
        else:
            min_max_recall[i] = (min_max_recall[i]-min_recall)/(max_recall-min_recall)
    
    normalised_recall[j][:] = min_max_recall
    general_recall[j][:] = recall
    heat_map_recall[j][:] = rank
    j = j+1

print(general_recall)

#Heatmap settings time
colors = ['#d6362b','#7373fa','#3535de','#07f727']
boundaries = [0,0.1,0.5,0.9,1]
cmap = LinearSegmentedColormap.from_list("custom", colors, N=len(colors))
black_out = LinearSegmentedColormap.from_list("custom", ['White','White','White'], N=3)
norm_black = BoundaryNorm([90,100,110], cmap.N, clip=True)
norm = BoundaryNorm(boundaries, cmap.N, clip=True)

ax = sns.heatmap(normalised_recall, cmap=cmap, norm=norm,linewidths=.5,linecolor='black',
                 cbar_kws = {})
sns.heatmap(normalised_recall,cmap=black_out,vmin = 99, vmax = 101, 
            mask=normalised_recall < 90, ax=ax,cbar=False,
            norm=norm_black,linecolor='black',linewidths=.5)
c_bar = ax.collections[0].colorbar
c_bar.set_ticks([0.05, 0.3, 0.7,0.95])
c_bar.set_ticklabels(['Bottom 10%', '50% - 90%', 
                      '10% - 50%','Top 10%'])
plt.title("Recall", fontsize=12)
plt.xlabel("Techniques", fontsize = 15)
ax.tick_params(axis='both', which='major')
ax.xaxis.set_ticks([x - y for x, y in zip(list(range(1,ranking_high+1)), ranking_high*[0.5])])
ax.set_xticklabels(method_names,rotation=90)
ax.yaxis.set_ticks([x - y for x, y in zip(list(range(1,len(onlyfiles)+1)), ranking_high*[0.5])])
ax.set_yticklabels(network_names,rotation=0, ha="right")
plt.ylabel("Data sets", fontsize = 15)
plt.tight_layout()
plt.show()


#---------------------------------Box plotting time----------------------------


test = pd.DataFrame({
    'Network type': network_names
})
print(test)



i = 0
for method in method_names:
    current_list = [k[i] for k in general_recall]
    # for j in range(len(current_list)):
    #     # if current_list[j] == -1:
    #     #     current_list[j] = 'NaN'
    test[method] = current_list
    i = i+1
print(test['EB_bootstrapping'])
test = pd.melt(test, ['Network type'])
#Set to -1 for precision
test = test[test.value != 100]


sns.set_theme(rc={'figure.figsize':(11.7,8.27)},style="whitegrid")
ax2 = sns.boxplot(x='variable', y='value', data = test,showfliers=False,
                 hue = 'variable'
                 ,fill=False,linewidth=1) 
ax2.set_xticklabels(test['variable'].unique(),rotation=90)
plt.xlabel("Techniques", fontsize = 15)
plt.ylabel("Recall", fontsize = 15)
ax2.grid(False)
plt.tight_layout()
plt.show()

