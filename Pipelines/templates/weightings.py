#This file will take in all the weightings in an attempt to get the average

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import stdev
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join

#Read in all the weightings
path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/weightings/All/Duds'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#Getting the method names
method_names = np.genfromtxt(path + '/' + onlyfiles[0], delimiter=',',dtype=None)[0]
print(onlyfiles)


#Setting up an empty dict
values = {method: [] for method in method_names}

#Now we put stuff in
for i in range(len(onlyfiles)):
    matrix = np.genfromtxt(path + '/' + onlyfiles[i], delimiter=',',dtype=None)
    for j in range(len(matrix[1])):
        values[matrix[0][j]].append(float(matrix[1][j]))

df = pd.DataFrame(columns = ['Type','-2std','Avg','+2std'])
i = 0
for key in values.keys():
    average = sum(values[key])/len(values[key])
    df.loc[i] = [key.replace(".csv", "").replace("GRN_","")
                 ,average-2*stdev(values[key]), average, average+2*stdev(values[key])]
    i = i+1

#Using this tutorial
#https://curbal.com/curbal-learning-portal/dumbbell-charts-in-matplotlib
fig, ax = plt.subplots(figsize=(10,6), facecolor = "white")

ax.grid(which="major", axis='both', color='#758D99', alpha=0.6, zorder=1)
ax.spines[['top','right','bottom']].set_visible(False)
ax.hlines(y=df['Type'], xmin=df['-2std'], xmax=df['+2std'], color='#758D99', zorder=2, linewidth=2, label='_nolegend_', alpha=.8)

ax.scatter(df['-2std'], df['Type'], label='-2std', s=60, color='#DB444B', zorder=3)
ax.scatter(df['Avg'], df['Type'], label='Average', s=60, color='#00a203', zorder=3)
ax.scatter(df['+2std'], df['Type'], label='+2std', s=60, color='#006BA2', zorder=3)

ax.legend(['-2std', 'Average', '+2std'], ncol=2, frameon=False, handletextpad=-.1, handleheight=1)
plt.yticks(fontsize=6)
plt.xlabel('Weighting of techniques')
plt.ylabel('Techniques')
plt.title('Dumbbell chart of the weights taken for Top 10 values')

plt.show()