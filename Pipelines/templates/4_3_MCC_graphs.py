#Because the 4.3 AUPR was too difficult to mofify, I will do the plotting here

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#This here olds the SERGIO multi files
path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/SERGIO/multi'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

name = []
AUPR_1 = []
AUPR_10 = []
AUPR_100 = []
AUPR = []

name = []
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
            if 'Matthew' in line:
                count = 1

    AUPR.append(AUPR_holder)


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
            if 'Matthew' in line:
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
            if 'Matthew' in line:
                count = 1
    
    AUPR_bee.append(AUPR_holder)

#--------------------------Now I can begin plotting everything------------------

for multi in [1,10,100]:
    #Holds the MCC values
    MCC_holder = []
    #Firstly get the SERGIO data
    i = 0
    for sergio in name:
        if sergio.endswith(str(multi)):
            MCC_holder = MCC_holder + AUPR[i][0:3]
        i = i+1
    
    #Now get the BEELINE stuff
    average_c = [0]*4
    average_m = [0]*4
    i = 0
    for x in name_bee:
        if x.endswith(str(multi)):
            if 'good' in x:
                average_c = [sum(j) for j in zip(average_c, AUPR[i])]
            if 'messy' in x:
                average_m = [sum(j) for j in zip(average_m, AUPR[i])]
            i = i+1
    
    bee_c_AUPR = [j/5 for j in average_c]
    bee_m_AUPR = [j/5 for j in average_m]
    MCC_holder = MCC_holder + bee_c_AUPR[0:3] + bee_m_AUPR[0:3]
    exec("MCC_" + str(multi) + " = MCC_holder")

#Get original data
average_orig_c = [0]*6
average_orig_m = [0]*6
MCC_non_bee = []
i = 0
for x in orig_name:
    if 'Clean test' in x:
        average_orig_c = [sum(j) for j in zip(average_orig_c, orig_AUPR[i])]
    elif 'Messy test' in x:
        average_orig_m = [sum(j) for j in zip(average_orig_m, orig_AUPR[i])]
    elif not x.endswith(str(100)):
        MCC_non_bee = MCC_non_bee + orig_AUPR[i][0:3]
    i = i+1
bee_orig_c_AUPR = [j/5 for j in average_orig_c][0:3]
bee_orig_m_AUPR = [j/5 for j in average_orig_m][0:3]

MCC_non_bee = MCC_non_bee + bee_orig_c_AUPR + bee_orig_m_AUPR

#Now I can put this into box plots

test = pd.DataFrame({
    'Inference method': ['No priors', 'Full priors', 'GENIE3 priors']*10,
    'Original': MCC_non_bee,
    'Multi = 1': MCC_1,
    'Multi = 10': MCC_10,
    'Multi = 100': MCC_100
})
test = pd.melt(test, ['Inference method'])
means = test.groupby("variable")["value"].mean()
means = [means[i] for i in [3,0,1,2]]
print(means)

sns.set_theme(rc={'figure.figsize':(11.7,8.27)},style="whitegrid")
ax = sns.boxplot(x='variable', y='value', data = test,showfliers=False,
                 hue = 'variable',fill=False,linewidth=2,
                 color='black',
                 showmeans=True,
                 meanprops={
                     'marker': 'o',
                     'markerfacecolor': 'red',
                     'markeredgecolor':'black'
                 }
                 ) 
ax.grid(False)
sns.stripplot(x='variable', y='value', data = test, ax=ax, hue = 'Inference method')
height = [-0.02]*4
for xtick in ax.get_xticks():
    ax.text(xtick,round(means[xtick],4) + height[xtick],round(means[xtick],4), 
            horizontalalignment='center',size='x-small',color='black',weight='semibold')
plt.legend(title="Inference methods")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel('multi',fontsize=15)
plt.ylabel("MCC score",fontsize=15)
# plt.title("MCC scores of various techniques using different 'multi' values",
#           fontsize = 20)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.show()







