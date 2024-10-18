#This plot here aims to get a nice dumbbell plot for section 4.4

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean

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
            if line == '\n' and count != 2:
                count = 0
            elif line == '\n' and count == 2:
                count = 1
            elif count == 1:
                AUPR_holder.append(float(re.findall("\d+\.\d+", line)[0]))
            if 'MCC' in line:
                count = 2
    
    AUPR.append(AUPR_holder)

#Get the SERGIO data
index_loc = []
labels = ['250 clean', '500 clean', '1000 clean', '2000 clean',
          '250 messy', '500 messy', '1000 messy', '2000 messy']
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

labels_dict = labels[0:4] + ['Beeline clean'] + labels[4:8] + ['Beeline messy']

#Now time to put it into a pandas dict
x = 10
type_dict = [type[0]]*x + [type[1]]*x + [type[2]]*x + [type[3]]*x + [type[4]]*x

MCC = []
for i in range(0,5):
    for j in range(0,8):
        if j == 4:
            MCC.append(bee_c_AUPR[i])
        MCC.append(AUPR[index_loc[j]][i])
    MCC.append(bee_m_AUPR[i])

min_val = []
max_val = []
avg = []
for i in range(0,5):
    min_val.append(min(MCC[10*i:(i+1)*10]))
    max_val.append(max(MCC[10*i:(i+1)*10]))
    avg.append(mean(MCC[10*i:(i+1)*10]))

test = pd.DataFrame({
    'Inference method': type,
    'Min': min_val,
    'Average': avg,
    'Max': max_val,
})

print(test)

fig, ax = plt.subplots(figsize=(10,6), facecolor = "white")
ax.grid(which="major", axis='both', color='#758D99', alpha=0.6, zorder=1)
ax.spines[['top','right','bottom']].set_visible(False)
ax.hlines(y=test['Inference method'], xmin=test['Min'], xmax=test['Max'], color='#758D99', zorder=2, linewidth=2, label='_nolegend_', alpha=.8)

ax.scatter(test['Min'], test['Inference method'], label='Min', s=60, color='#DB444B', zorder=3)
ax.scatter(test['Max'], test['Inference method'], label='Max', s=60, color='#DB444B', zorder=3)

plt.xlabel('MCC score')
plt.ylabel('Techniques')


plt.show()




