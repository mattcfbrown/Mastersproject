#This creates the figure I wish to use for section 4.2

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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
    'Dataset': [labels[x] for x in [1,2,3,4]] + ['Beeline Clean Average'] +  
                   [labels[x] for x in [6,7,8,9]] +['Beeline Messy Average'],
    'No priors': [[a[0] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4]] +
        [bee_c_AUPR[0]] + [[a[0] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [6,7,8,9]]
        + [bee_m_AUPR[0]], 
    'Full priors': [[a[1] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4]] +
        [bee_c_AUPR[1]] + [[a[1] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [6,7,8,9]]
        + [bee_m_AUPR[1]], 
    'GENIE3 priors': [[a[2] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4]] +
        [bee_c_AUPR[2]] + [[a[2] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [6,7,8,9]]
        + [bee_m_AUPR[2]], 
    # 'NLNET priors': [[a[3] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4,6,7,8,9]] +
    #     [bee_c_AUPR[3], bee_m_AUPR[3]],
    'GENIE3': [[a[4] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4]] +
        [bee_c_AUPR[4]] + [[a[4] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [6,7,8,9]]
        + [bee_m_AUPR[4]], 
    'PUC': [[a[5] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [1,2,3,4]] +
        [bee_c_AUPR[5]] + [[a[5] for a in list(map(AUPR.__getitem__, index_loc))][x] for x in [6,7,8,9]]
        + [bee_m_AUPR[5]], 
})
test = pd.melt(test, ['Dataset'])
means = test.groupby("variable")["value"].mean()
means = [means[i] for i in [3,0,2,1,4]]

rng = np.random.default_rng()
jitter = 0.1*rng.normal(size=50)
print(jitter)
values = [0]*10+ [1]*10+ [2]*10+ [3]*10+ [4]*10
type = ['clean']*5 + ['messy'] *5
type = type * 5
values = [sum(i) for i in zip(values,jitter)]
test.insert(3,"X_loc",values)
test.insert(4,'Type',type)

# print(means)

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
ax = sns.boxplot(x='variable', y='value', data = test,showfliers=False
                 ,fill=False,linewidth=1,
                 color = 'black',
                 showmeans=True,
                 meanprops={
                     'marker': 'o',
                     'markerfacecolor': 'red',
                     'markeredgecolor':'black'
                 }) 
ax.grid(False)

# clean = pd.DataFrame(columns=['Data_type', 'variable', 'value'])
# for type in labels[0:5] + ['Beeline Clean Average']:
#     clean = pd.concat([clean ,test[test['Data_type'] == type]],ignore_index=True)

# messy = pd.DataFrame(columns=['Data_type', 'variable', 'value'])
# for type in labels[5:10] + ['Beeline Messy Average']:
#     messy = pd.concat([messy ,test[test['Data_type'] == type]],ignore_index=True)

# i = 0
# markers = ["s","p","P","*","H"]
# for type in labels[1:5] + ['Beeline Clean Average']:
#     print(type)
#     clean = pd.DataFrame(columns=['Data_type', 'variable', 'value'])
#     clean = pd.concat([clean ,test[test['Data_type'] == type]],ignore_index=True)
#     sns.swarmplot(x='variable', y='value', data = clean, ax=ax, hue = 'Data_type',
#                 palette = ['green'],marker = str(markers[i]),
#                 size = 7, dodge = True)
#     i = i+1
# i = 0
# for type in labels[6:10] + ['Beeline Messy Average']:
#     print(type)
#     clean = pd.DataFrame(columns=['Data_type', 'variable', 'value'])
#     clean = pd.concat([clean ,test[test['Data_type'] == type]],ignore_index=True)
#     sns.swarmplot(x='variable', y='value', data = clean, ax=ax, hue = 'Data_type',
#                 palette = ['red'],marker = str(markers[i]),
#                 size = 7, dodge = True)
    # i = i+1
                # palette = sns.color_palette(palette='Greens_d', n_colors=5) + 
                #     sns.color_palette(palette='Reds', n_colors=5))

# markers = ["s","p","P","*","H"]
# print(test[test['Data_type'] == 'Clean 250'])
# sns.scatterplot(x='variable', y='value', data = test[test['Data_type'] == 'Clean 250'],
#                ax=ax, hue = 'Data_type',
#                 palette = ['green'],marker = markers,
#                 size = 7)


sns.scatterplot(data = test,
                x = 'X_loc', 
                y = 'value',
                hue = (['Clean']*5 + ['Messy']*5)*5,
                style = ['250 cell', '500 cell', '1000 cell', '2000 cell',
                         'Beeline']*10,
                ax=ax,
                palette=['#28d73c','#d6293b'],
                s=100)
# i = 0
# markers = ["s","p","P","*","H"]
# for type in labels[1:5] + ['Beeline Clean Average']:
#     clean = pd.DataFrame(columns=['Data_type', 'variable', 'value', 'X_loc'])
#     clean = pd.concat([clean ,test[test['Data_type'] == type]],ignore_index=True)
#     sns.scatterplot(data = clean,
#                     x = 'X_loc', 
#                     y = 'value',
#                     marker = markers[0],
#                     palette = ["green"],
#                     ax=ax)
#     i = i+1
    


height = [-0.02]*4 + [0.01]
for xtick in ax.get_xticks():
    ax.text(xtick,round(means[xtick],4) + height[xtick],round(means[xtick],4), 
            horizontalalignment='center',size='x-small',color='black',weight='semibold')
plt.legend(title="Datasets"
           )
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.xlabel("Network Inference techniques",fontsize=15)
plt.ylabel("MCC score",fontsize=15)
# plt.title("MCC scores of various techniques using different datasets",
#           fontsize = 20)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.show()
# print('\n\n')

print(test)