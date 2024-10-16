#Graph which is able to get the data I need
#Can do this tomorrow :)

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/SERGIO/multi'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

name = []
AUPR_1 = []
AUPR_10 = []
AUPR_100 = []
AUPR = []
for file in onlyfiles:
    name_hold = file.replace("matthew_", "")
    name_hold = name_hold.replace(".txt", "")
    name_hold = name_hold.replace("_Test_0", "")
    name_hold = name_hold.replace("_", " ")
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
    # if ' 100' in name_hold:
    #     AUPR_100.append(AUPR_holder)
    # elif ' 10' in name_hold:
    #     AUPR_10.append(AUPR_holder)
    # else:
    #     AUPR_1.append(AUPR_holder)

#Now get the original stuff
#----------------------Original---------------------------------------------
orig_path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/results/beeline/metrics'
onlyfiles = [f for f in listdir(orig_path) if isfile(join(orig_path, f))]


orig_name = []
orig_AUPR = []
for file in onlyfiles:
    name_hold = file.replace("matthew_25_", "")
    name_hold = name_hold.replace("matthew_11_", "")
    name_hold = name_hold.replace(".txt", "")
    name_hold = name_hold.replace("_Test_0", "")
    if 'Messy' not in name_hold:
        name_hold = 'Clean_' + name_hold
    name_hold = name_hold.replace("_", " ")
    orig_name.append(name_hold)

    AUPR_holder = []
    count = 0
    #Now extract everything
    with open(orig_path + '/' + file,'r') as f:
        for line in f:
            if line == '\n':
                count = 0
            if count == 1:
                AUPR_holder.append(float(re.findall("\d+\.\d+", line)[0]))
            if 'AUPR' in line:
                count = 1
    
    orig_AUPR.append(AUPR_holder)

#--------------------------------Now get the BEELINE stuff--------------------------
beeline_path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/BEELINE/multi'
beefiles = [f for f in listdir(beeline_path) if isfile(join(beeline_path, f))]
name_bee = []
AUPR_bee = []
for file in beefiles:
    name_hold = file.replace("matthew_", "")
    name_hold = name_hold.replace(".txt", "")
    name_hold = name_hold.replace("_Test_0", "")
    name_hold = name_hold.replace("_", " ")
    name_bee.append(name_hold)

    AUPR_holder = []
    count = 0
    #Now extract everything
    with open(beeline_path + '/' + file,'r') as f:
        for line in f:
            if line == '\n':
                count = 0
            if count == 1:
                AUPR_holder.append(float(re.findall("\d+\.\d+", line)[0]))
            if 'AUPR' in line:
                count = 1
    
    AUPR_bee.append(AUPR_holder)

#--------------------------Now I can begin plotting everything------------------
#Get original BEELINE
average_orig_c = [0]*6
average_orig_m = [0]*6
i = 0
for x in orig_name:
    if 'Clean test' in x:
        average_orig_c = [sum(j) for j in zip(average_orig_c, orig_AUPR[i])]
    if 'Messy test' in x:
        average_orig_m = [sum(j) for j in zip(average_orig_m, orig_AUPR[i])]
    i = i+1
bee_orig_c_AUPR = [j/5 for j in average_orig_c]
bee_orig_m_AUPR = [j/5 for j in average_orig_m]


labels = ['Clean 250', 'Clean 500', 'Clean 1000', 'Clean 2000',
          'Messy 250', 'Messy 500', 'Messy 1000', 'Messy 2000']


for multi in [1,10,100]:

    average_c = [0]*4
    average_m = [0]*4
    i = 0
    #Step 1 Get BEELINE averages
    for x in name_bee:
        if x.endswith(str(multi)):
            if 'good' in x:
                average_c = [sum(j) for j in zip(average_c, AUPR[i])]
            if 'messy' in x:
                average_m = [sum(j) for j in zip(average_m, AUPR[i])]
            i = i+1
        
    bee_c_AUPR = [j/5 for j in average_c]
    bee_m_AUPR = [j/5 for j in average_m]

    #Now we get their subtraction
    bee_c_diff = [a-b for a,b in zip(bee_c_AUPR[0:3],bee_orig_c_AUPR[0:3])]
    bee_m_diff = [a-b for a,b in zip(bee_m_AUPR[0:3],bee_orig_m_AUPR[0:3])]

    #Now we get the differences
    diff = []
    for lab in labels:
        #Get the label
        if 'Clean' in lab:
            value = str(re.findall("\d+", lab)[0]) + 'cells ' + str(multi)
        else:
            value = str(re.findall("\d+", lab)[0]) + 'cells messy ' + str(multi)
        #Find it's location
        orig_location = orig_name.index(lab)
        location = name.index(value)
        # print(location)
        #Difference
        differ = [a-b for a,b in zip(AUPR[location][0:3],orig_AUPR[orig_location][0:3])]
        diff.append(differ)
        
    #Now we can put this data into the pandas database
    type = ['No priors', 'Full priors', 'GENIE3 priors']
    d = []
    for p in range(8):
        for j in range(3):
            d.append(
                {
                    'Name': labels[p],
                    'Difference': diff[p][j],
                    'technique':  type[j]
                }
            )
    
    #BEELINE time
    for j in range(3):
        d.append(
            {
                'Name': 'Beeline Clean Average',
                'Difference': bee_c_diff[j],
                'technique':  type[j]
            }
        )

    for j in range(3):
        d.append(
            {
                'Name': 'Beeline Messy Average',
                'Difference': bee_m_diff[j],
                'technique':  type[j]
            }
        )
    
    data = pd.DataFrame(d)
    print(data)
    
    #It's plotting time
    sns.set_theme(style="darkgrid")
    sns.barplot(data=data, x='Name', y='Difference', hue='technique')
    title = 'Difference of AUPR between normal prior and alternate prior with multi = ' + str(multi)
    plt.title(title)
    plt.xlabel("Datasets")
    plt.xticks(rotation=30,fontsize=6,horizontalalignment='right')
    plt.show()
        



