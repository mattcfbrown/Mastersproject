#This creates the figure I wish to use for section 4.3

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from statistics import mean

path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/SERGIO'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]


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

#We get the location of the two
loc = ['250cells', '500cells', '1000cells', '2000cells']
index_loc = []
difference = []
#For all 4 we are getting
    #clean 2.2
    #clean 2.944
    #messy 2.2
    #messy 2.944
for item in loc:
    index_loc.append(name.index(item + ' 2.2'))
    index_loc.append(name.index(item + ' 2.944'))
    index_loc.append(name.index(item + ' messy 2.2'))
    index_loc.append(name.index(item + ' messy 2.944'))

difference = []
for i in range(8):
    difference.append([a-b for a,b in zip(AUPR[index_loc[2*i]],
                              AUPR[index_loc[2*i + 1]])])

type = ['No priors', 'Full priors', 'GENIE3 priors']
dataset = ['250cells', '250cells messy',
           '500cells', '500cells messy',
           '1000cells', '1000cells messy',
           '2000cells', '2000cells messy']

d = []
for p in range(8):
    for j in range(3):
        d.append(
            {
                'Name': dataset[p],
                'Difference': difference[p][j],
                'technique':  type[j]
            }
        )

#-------Beeline stuff-----
beeline_path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/BEELINE'
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

    

average_c_22 = [0]*3
average_c_2944 = [0]*3
average_m_22 = [0]*3
average_m_2944 = [0]*3
i = 0
for x in name_bee:
    if 'good' in x:
        if '2.2' in x:
            average_c_22 = [sum(j) for j in zip(average_c_22, AUPR_bee[i])]
        else:
            average_c_2944 = [sum(j) for j in zip(average_c_2944, AUPR_bee[i])]
    if 'messy' in x:
        if '2.2' in x:
            average_m_22 = [sum(j) for j in zip(average_m_22, AUPR_bee[i])]
        else:
            average_m_2944 = [sum(j) for j in zip(average_m_2944, AUPR_bee[i])]
    i = i+1
bee_c_AUPR = [ a-b for a,b in zip([j/5 for j in average_c_22],
                                  [j/5 for j in average_c_2944])]
bee_m_AUPR = [ a-b for a,b in zip([j/5 for j in average_m_22],
                                  [j/5 for j in average_m_2944])]


for j in range(3):
    d.append(
        {
            'Name': 'Beeline Clean Average',
            'Difference': bee_c_AUPR[j],
            'technique':  type[j]
        }
    )

for j in range(3):
    d.append(
        {
            'Name': 'Beeline Messy Average',
            'Difference': bee_m_AUPR[j],
            'technique':  type[j]
        }
    )



#-------------------------



test = pd.DataFrame({
    'Inference method': type*10,
    '2.2': list(itertools.chain(
            *[j[0:3] for j in
            [AUPR[index_loc[2*i]] for i in range(8)]
            ])) + [k/5 for k in average_c_22] + [k/5 for k in average_m_22],
    '2.944': list(itertools.chain(
            *[j[0:3] for j in
            [AUPR[index_loc[2*i + 1]] for i in range(8)]
            ])) + [k/5 for k in average_c_2944] + [k/5 for k in average_m_2944]
})
test = pd.melt(test, ['Inference method'])
means = test.groupby("variable")["value"].mean()
print(means)


# data = pd.DataFrame(d)
# print(data)


# sns.set_theme(style="darkgrid")
# sns.barplot(data=data, x='Name', y='Difference', hue='technique')
# plt.title('AUPR difference between 2.2 and 2.944 of various techniques')
# plt.xlabel("Datasets")
# plt.xticks(rotation=30,fontsize=6,horizontalalignment='right')
# plt.show()

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
sns.stripplot(x='variable', y='value', data = test, ax=ax, hue = 'Inference method')
for xtick in ax.get_xticks():
    ax.text(xtick,round(means[xtick],4) - 0.01,round(means[xtick],4), 
            horizontalalignment='center',size='x-small',color='red',weight='semibold')
plt.legend(title="Inference methods")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel(r'$\omega^{0}$',fontsize=15)
plt.ylabel("Matthew's correlation coefficient scores",fontsize=15)
plt.title("MCC scores of various techniques using different " + r'$\omega^{0}$' + " values",
          fontsize = 20)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.show()
