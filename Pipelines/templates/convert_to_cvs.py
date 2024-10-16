#This here just takes all the files and puts them into something readable

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join
from statistics import mean

path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/BEELINE'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#We want to split them up
beeline = []
for file in onlyfiles:
    if '10' in file and ('good' in file or 'messy' in file) and '100' not in file:
        beeline.append(file)

fi_no_clean = []
mcc_no_clean  = []
aupr_no_clean  = []
fi_full_clean  = []
mcc_full_clean  = []
aupr_full_clean  = []
fi_no_messy = []
mcc_no_messy = []
aupr_no_messy  = []
fi_full_messy = []
mcc_full_messy = []
aupr_full_messy  = []
section = []
for file in beeline:
    section.append(file)
    if 'messy' in file:
        clean = 0
    else:
        clean = 1
    with open(path + '/' + file,'r') as f:
        count = 0
        counter = 0
        for line in f:
            if 'No priors' in line and line != 'No priors:\n':
                if count == 0:
                    if clean == 1:
                        mcc_no_clean.append(float(re.findall("\d+\.\d+", line)[0]))
                    else:
                        mcc_no_messy.append(float(re.findall("\d+\.\d+", line)[0]))
                    count = 1
                elif count == 1:
                    if clean == 1:
                        fi_no_clean.append(float(re.findall("\d+\.\d+", line)[0]))
                    else:
                        fi_no_messy.append(float(re.findall("\d+\.\d+", line)[0]))
                    count = 2
                else:
                    if clean == 1:
                        aupr_no_clean.append(float(re.findall("\d+\.\d+", line)[0]))
                    else:
                        aupr_no_messy.append(float(re.findall("\d+\.\d+", line)[0]))
                    count = 0
            if 'Full priors' in line and line != 'Full priors:\n':
                if counter == 0:
                    if clean == 1:
                        mcc_full_clean.append(float(re.findall("\d+\.\d+", line)[0]))
                    else:
                        mcc_full_messy.append(float(re.findall("\d+\.\d+", line)[0]))
                    counter = 1
                elif counter == 1:
                    if clean == 1:
                        fi_full_clean.append(float(re.findall("\d+\.\d+", line)[0]))
                    else:
                        fi_full_messy.append(float(re.findall("\d+\.\d+", line)[0]))
                    counter = 2
                else:
                    if clean == 1:
                        aupr_full_clean.append(float(re.findall("\d+\.\d+", line)[0]))
                    else:
                        aupr_full_messy.append(float(re.findall("\d+\.\d+", line)[0]))
                    counter = 0



print(mean(mcc_no_clean))
print(mean(mcc_full_clean))
print(mean(fi_no_clean))
print(mean(fi_full_clean))
print(mean(aupr_no_clean))
print(mean(aupr_full_clean))
print(mean(mcc_no_messy))
print(mean(mcc_full_messy))
print(mean(fi_no_messy))
print(mean(fi_full_messy))
print(mean(aupr_no_messy))
print(mean(aupr_full_messy))



