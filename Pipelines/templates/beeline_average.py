#Just gonna get averages of beeline EB outputs

import re
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join

# path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/PIDC_testing/results/thesis/metrics/prior'
# path ='/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/results/beeline/metrics'
# path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/BEELINE'
path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/metrics/BEELINE/multi'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

average_AUPR_c = [0]*6
average_MCC_c  = [0]*6
average_F1_c   = [0]*6
average_AUPR_m = [0]*6
average_MCC_m  = [0]*6
average_F1_m  = [0]*6
average_found       = [0]*6
average_missed      = [0]*6
average_incorrect   = [0]*6

matthew = 0
f1 = 0
AUPR = 0
for file in onlyfiles:
    found_holder = []
    missed_holder = []
    incorrect_holder = []
    with open(path + '/' + file,'r') as f:
        if '100.txt' in file:
            for line in f:

                #Sets up the data
                if line == '\n':
                    if (matthew == 1):
                        if 'messy' not in file:
                            average_MCC_c = [sum(i) for i in zip(average_MCC_c, holder)]
                        else:
                            average_MCC_m = [sum(i) for i in zip(average_MCC_m, holder)]
                    if (f1 == 1):
                        if 'messy' not in file:
                            average_F1_c = [sum(i) for i in zip(average_F1_c, holder)]
                        else:
                            average_F1_m = [sum(i) for i in zip(average_F1_m, holder)]
                    if (AUPR == 1):
                        if 'messy' not in file:
                            average_AUPR_c = [sum(i) for i in zip(average_AUPR_c, holder)]
                        else:
                            average_AUPR_m = [sum(i) for i in zip(average_AUPR_m, holder)]
                    matthew = 0
                    f1 = 0
                    AUPR = 0
                
                #Puts in the values
                if matthew == 1:
                    holder.append(float(re.findall("\d+\.\d+", line)[0]))
                if f1 == 1:
                    holder.append(float(re.findall("\d+\.\d+", line)[0]))
                if AUPR == 1:
                    holder.append(float(re.findall("\d+\.\d+", line)[0]))


                #Get what the next line should be
                if 'Matthew' in line:
                    holder = []
                    matthew = 1
                if 'F1' in line:
                    holder = []
                    f1 = 1
                if 'AUPR' in line:
                    holder = []
                    AUPR = 1
                
                if 'Agreed upon gene connections' in line:
                    found_holder.append(float(re.findall("\d+", line)[0]))
                if 'Missed gene connections' in line:
                    missed_holder.append(float(re.findall("\d+", line)[0]))
                if 'Incorrect gene connections' in line:
                    average_incorrect.append(float(re.findall("\d+", line)[0]))
    
    average_found = [sum(i) for i in zip(average_found, found_holder)]
    average_missed = [sum(i) for i in zip(average_found, missed_holder)]
    average_incorrect = [sum(i) for i in zip(average_found, average_incorrect)]
    

print('Stats info \n')
# print([i/5 for i in average_MCC_c])
# print([i/5 for i in average_MCC_m])
print([i/5 for i in average_F1_c])
print('\n')
print([i/5 for i in average_F1_m])
print('\n')
print([i/5 for i in average_AUPR_c])
print('\n')
print([i/5 for i in average_AUPR_m])
# print('\nEdge info \n')
# print([i/5 for i in average_found])
# print([i/5 for i in average_missed])
# print([i/5 for i in average_incorrect])
            