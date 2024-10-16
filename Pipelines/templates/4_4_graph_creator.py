#This here is the code needed for chapter 4.4

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/PIDC_testing/results/thesis/metrics/prior'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

name = []
AUPR = []
for file in onlyfiles:
    name_hold = file.replace("metrics_", "")
    name_hold = name_hold.replace(".txt", "")
    name_hold = name_hold.replace("_", " ")
    name_hold = name_hold.replace("cells", "")
    name_hold = name_hold.replace("good", "")
    if 'messy' not in name_hold:
        name_hold = name_hold + ' clean'
    name.append(name_hold)

    AUPR_holder = []
    count = 0
    #Now extract everything
    with open(path + '/' + file,'r') as f:
        for line in f:
            if line == '\n':
                count = 0
            if count == 1:
                AUPR_holder.append(float(re.findall("\d+\.\d+", line)[0]))
            if 'AUPR' in line:
                count = 1
    
    AUPR.append(AUPR_holder)

#Get the SERGIO data
index_loc = []
labels = ['100 clean', '250 clean', '500 clean', '1000 clean', '2000 clean',
          '100 messy', '250 messy', '500 messy', '1000 messy', '2000 messy']
for x in labels:
    index_loc.append(name.index(x))

type = ['No priors', 'Full priors', 'GENIE3 priors', 'GENIE3', 'PIDC']

average_c = [0]*5
average_m = [0]*5
i = 0
for x in name:
    if 'clean' in x and 'bee' in x:
        average_c = [sum(j) for j in zip(average_c, AUPR[i])]
    if ('messy' in x) and ('bee' in x):
        average_m = [sum(j) for j in zip(average_m, AUPR[i])]
    i = i+1
bee_c_AUPR = [j/5 for j in average_c]
bee_m_AUPR = [j/5 for j in average_m]

d = []
for p in index_loc:
    for j in range(5):
        d.append(
            {
                'Name': name[p],
                'AUPR': AUPR[p][j],
                'technique':  type[j]
            }
        )

for j in range(5):
    d.append(
        {
            'Name': 'Beeline Clean Average',
            'AUPR': bee_c_AUPR[j],
            'technique':  type[j]
        }
    )

for j in range(5):
    d.append(
        {
            'Name': 'Beeline Messy Average',
            'AUPR': bee_m_AUPR[j],
            'technique':  type[j]
        }
    )

data = pd.DataFrame(d)

#Now graphing time
sns.set_theme(style="darkgrid")
sns.barplot(data=data, x='Name', y='AUPR', hue='technique')
plt.title('AUPR scores of various techniques using PIDC test statisitcs and altered prior')
plt.xlabel("Datasets")
plt.xticks(rotation=30,fontsize=6,horizontalalignment='right')
plt.show()



