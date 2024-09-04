#This here will aim to generate three graphs for use
    #1: Heat map
    #2: Bar graph
    #3: Table of data

#We firstly need the following libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import seaborn as sns
from statistics import mean
import pandas as pd
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join

#Firstly read in all files
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


#--------------------------------PART 1: GETTING THE HEATMAP-----------------------------------
#We need to create a ranking for each method
#Gets the ranking number
ranking_high = len(method_names)
heat_map = np.zeros((len(onlyfiles),
                     ranking_high)) #Will hold the heat_map values 
j = 0
for row in AUPR:
    ranks = list(np.argsort((-1*row).argsort()))
    #Editing the rank value
    loc_before = ranks.index(0)
    for i in range(1,ranking_high):
        loc = ranks.index(i)
        #If we have the same value, update the rankings
        if row[loc] == row[loc_before]:
            ranks[loc] = ranks[loc_before]
        loc_before = loc
    
    heat_map[j][:] = ranks
    j = j+1

#Heatmap settings time
colors = ['#07f727','#7373fa','#3535de','#d6362b']
boundaries = [0,1,7,20,29]
cmap = LinearSegmentedColormap.from_list("custom", colors, N=len(colors))
norm = BoundaryNorm(boundaries, cmap.N, clip=True)

ax = sns.heatmap(heat_map, cmap=cmap, norm=norm,linewidths=.5,linecolor='black')
plt.title("AUPR score heatmap", fontsize=12)
plt.xlabel("Techniques", fontsize=8)
ax.tick_params(axis='both', which='major', labelsize=5)
ax.xaxis.set_ticks([x - y for x, y in zip(list(range(1,ranking_high+1)), ranking_high*[0.5])])
ax.set_xticklabels(method_names,rotation=45, ha="right")
ax.yaxis.set_ticks([x - y for x, y in zip(list(range(1,len(onlyfiles)+1)), ranking_high*[0.5])])
ax.set_yticklabels(network_names,rotation=0, ha="right")
plt.ylabel("Data sets", fontsize=8)
plt.show()

#--------------------------------PART 2: GETTING THE BARCHART-----------------------------------

#We want a list of average AUPR values
average_value = []
#Gets average value
for column in AUPR.T:
   average_value.append(mean(column))
#Gets sorted index locations
sorted_average = [i[0] for i in sorted(enumerate(average_value), key=lambda x:x[1])]

#Plot
sns.set_style('darkgrid')
ax_bar = sns.barplot(y=[method_names[i] for i in sorted_average], 
            x=sorted(average_value))
ax_bar.tick_params(axis='both', which='major', labelsize=7)
ax_bar.bar_label(ax_bar.containers[0], size = 8)
plt.title("Average AUPR scores for (n=17) networks", fontsize=15)
plt.xlabel("AUPR value", fontsize=10)
plt.ylabel("Methods", fontsize=10)
plt.show()

#--------------------------------PART 3: GETTING THE TABLE-----------------------------------

#We firstly create the dataframe

data = {
    'Method': [],
    'Average Precision': [],
    'Average Recall': [],
}
df = pd.DataFrame(data)

#We loop over it all
i = 0
for method in method_names:
    TP = sum(found_array[:,i])
    FN = sum(missed_array[:,i])
    FP = sum(incorrect_array[:,i])
    new_row = {'Method': method, 'Average Precision': TP/(TP+FP), 'Average Recall': TP/(TP+FN)}
    df.loc[len(df)] = new_row
    i = i+1

print(df)
