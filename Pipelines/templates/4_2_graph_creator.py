#This creates the figure I wish to use for section 4.2

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/results/beeline/metrics'
# path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/results/SERGIO/metrics'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]


name = []
AUPR = []
for file in onlyfiles:
    name_hold = file.replace("matthew_25_", "")
    name_hold = name_hold.replace("matthew_11_", "")
    name_hold = name_hold.replace(".txt", "")
    name_hold = name_hold.replace("_Test_0", "")
    if 'Messy' not in name_hold:
        name_hold = 'Clean_' + name_hold
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


#Get the SERGIO data
index_loc = []
labels = ['Clean 100', 'Clean 250', 'Clean 500', 'Clean 1000', 'Clean 2000',
          'Messy 100', 'Messy 250', 'Messy 500', 'Messy 1000', 'Messy 2000']
for x in labels:
    index_loc.append(name.index(x))

type = ['No priors', 'Full priors', 'GENIE3 priors', 'NLNET priors', 'GENIE3', 'PUC']

average_c = [0]*6
average_m = [0]*6
i = 0
for x in name:
    if 'Clean test' in x:
        average_c = [sum(j) for j in zip(average_c, AUPR[i])]
    if 'Messy test' in x:
        average_m = [sum(j) for j in zip(average_m, AUPR[i])]
    i = i+1
bee_c_AUPR = [j/5 for j in average_c]
bee_m_AUPR = [j/5 for j in average_m]




d = []
for p in index_loc:
    for j in range(6):
        d.append(
            {
                'Name': name[p],
                'AUPR': AUPR[p][j],
                'technique':  type[j]
            }
        )

for j in range(6):
    d.append(
        {
            'Name': 'Beeline Clean Average',
            'AUPR': bee_c_AUPR[j],
            'technique':  type[j]
        }
    )

for j in range(6):
    d.append(
        {
            'Name': 'Beeline Messy Average',
            'AUPR': bee_m_AUPR[j],
            'technique':  type[j]
        }
    )



data = pd.DataFrame(d)
test = pd.DataFrame({
    'Data_type': [labels[x] for x in [1,2,3,4,6,7,8,9]] + ['Beeline Clean Average', 'Beeline Messy Average'],
    'No priors': [[a[0] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4,6,7,8,9]] +
        [bee_c_AUPR[0], bee_m_AUPR[0]], 
    'Full priors': [[a[1] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4,6,7,8,9]] +
        [bee_c_AUPR[1], bee_m_AUPR[1]],
    'GENIE3 priors': [[a[2] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4,6,7,8,9]] +
        [bee_c_AUPR[2], bee_m_AUPR[2]],
    'NLNET priors': [[a[3] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4,6,7,8,9]] +
        [bee_c_AUPR[3], bee_m_AUPR[3]],
    'GENIE3': [[a[4] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4,6,7,8,9]] +
        [bee_c_AUPR[4], bee_m_AUPR[4]],
    'PUC': [[a[5] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4,6,7,8,9]] +
        [bee_c_AUPR[5], bee_m_AUPR[5]]
})


test = pd.melt(test, ['Data_type'])
means = test.groupby("variable")["value"].mean()
means = [means[i] for i in [4,0,2,3,1,5]]
print(means)

# print(type)
# print(
#     [a[0] for a in 
#      list(map(AUPR.__getitem__, index_loc))
#      ]
#     )
# print(AUPR[index_loc[0]])
# print(test)

#Now graphing time
# sns.set_theme(style="whitegrid")
# sns.barplot(data=data, x='Name', y='AUPR', hue='technique')
# plt.title('AUPR scores of various techniques')
# plt.xlabel("Datasets")
# plt.xticks(rotation=30,fontsize=6,horizontalalignment='right')
# plt.show()

# sns.lineplot(data=test, x='Data_type', y='value', hue='variable')
# plt.title('AUPR scores of various techniques')
# plt.xlabel("Datasets")
# plt.xticks(rotation=30,fontsize=6,horizontalalignment='right')
# plt.show()

# sns.violinplot(data=test, x='variable', y='value')
sns.set_theme(rc={'figure.figsize':(11.7,8.27)},style="whitegrid")
ax = sns.boxplot(x='variable', y='value', data = test,showfliers=False,
                 hue = 'variable',fill=False,linewidth=2,
                 showmeans=True,
                 meanprops={
                     'marker': 'o',
                     'markerfacecolor': 'red',
                     'markeredgecolor':'black'
                 }) 
ax.grid(False)
sns.stripplot(x='variable', y='value', data = test, ax=ax, hue = 'Data_type')
height = [-0.02]*5 + [0.01]
for xtick in ax.get_xticks():
    ax.text(xtick,round(means[xtick],4) + height[xtick],round(means[xtick],4), 
            horizontalalignment='center',size='x-small',color='black',weight='semibold')
plt.legend(title="Datasets")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel("Network Inference techniques",fontsize=15)
plt.ylabel("Matthew's correlation coefficient scores",fontsize=15)
plt.title("MCC scores of various techniques using different datasets",
          fontsize = 20)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.show()