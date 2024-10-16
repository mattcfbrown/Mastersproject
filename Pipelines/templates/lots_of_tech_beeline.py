#Hopefully will be the last script I have to write

#This here takes in the beeline reads and returns a consensus text file.

#Should be the Average AUPR AND total connections found


#All the libraries needed
import numpy as np
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join


#Firstly read in all files
path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/metrics/Beeline'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#We get the names of the methods we wish to test:
method_names = []
file = open(path + '/' + onlyfiles[1], 'r')
for line in file.readlines():
    if 'connection' in line:
        break
    if len(line.rstrip()) > 0 and 'GRNVBEM' not in line and 'AUPR' not in line:
       method_names.append(line.split(':', 1)[0])

#We now want to get the data from each file
network_names = []                       #Holds the names of the networks
AUPR = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the AUPR scores
found_array = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the found scores
missed_array = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the missed scores
incorrect_array = np.zeros((len(onlyfiles),
                     len(method_names))) #Will hold the incorrect scores
i = 0

for method in onlyfiles:
    #We firstly get all the names
    name = method.replace("Metrics_", "")
    name = name.replace(".txt", "")
    network_names.append(name)
    
    #We are now getting all the values needed
    j = 0
    AUPR_done = 0
    found = 0
    missed = 0
    incorrect = 0
    file = open(path + '/' + method, 'r')
    for line in file.readlines():
        if 'connection' in line:
            AUPR_done = 1
        if AUPR_done == 0:
            if len(line.rstrip()) > 0 and 'GRNVBEM' not in line and 'AUPR' not in line:
                value = line.replace(method_names[j] + ': ', '')
                value = value.replace('\n', '')
                AUPR[i][j] = value
                j = j+1
        else:
            #Get the method
            if 'connection' in line and 'GRNVBEM' not in line:
                grnvbem = 0
                current_method = line.replace(' connection information\n', '')
                j = method_names.index(current_method)
            if 'GRNVBEM' in line:
                incorrect = 0
                found = 0
                missed = 0
                grnvbem = 1
            #Defines where we are
            if 'FOUND' in line and grnvbem == 0:
                found = 1
                incorrect = 0
            elif 'MISSED' in line and grnvbem == 0:
                missed = 1
                found = 0
            elif 'INCORRECT' in line and grnvbem == 0:
                incorrect = 1
                missed = 0
            #We now do the counting
            if found == 1 and '---' in line:
                found_array[i][j] = found_array[i][j] + 1
            elif missed == 1 and '---' in line:
                missed_array[i][j] = missed_array[i][j] + 1
            elif incorrect == 1 and '---' in line:
                incorrect_array[i][j] = incorrect_array[i][j] + 1
    
    i = i+1

#Now time to put this back into a text file
c_index = []
m_index = []
i = 0
for data in network_names:
    if 'good' in data:
        c_index.append(i)
    else:
        m_index.append(i)
    i = i+1

#Now it's text file time

holder = [0]*len(method_names)
f_holder = [0]*len(method_names)
m_holder = [0]*len(method_names)
i_holder = [0]*len(method_names)
for i in c_index:
    #Average out AUPR scores
    holder = [a+b for a,b in zip(holder,AUPR[i])]
    f_holder = [a+b for a,b in zip(f_holder,found_array[i])]
    m_holder = [a+b for a,b in zip(m_holder,missed_array[i])]
    i_holder = [a+b for a,b in zip(i_holder,incorrect_array[i])]
c_AUPR = [j/5 for j in holder]
c_found = f_holder
c_missed = m_holder
c_incorrect = i_holder


holder = [0]*len(method_names)
f_holder = [0]*len(method_names)
m_holder = [0]*len(method_names)
i_holder = [0]*len(method_names)
for i in m_index:
    #Average out AUPR scores
    holder = [a+b for a,b in zip(holder,AUPR[i])]
    f_holder = [a+b for a,b in zip(f_holder,found_array[i])]
    m_holder = [a+b for a,b in zip(m_holder,missed_array[i])]
    i_holder = [a+b for a,b in zip(i_holder,incorrect_array[i])]
m_AUPR = [j/5 for j in holder]
m_found = f_holder
m_missed = m_holder
m_incorrect = i_holder



with open('Metrics_beeline_clean.txt', 'w') as f:
        #Firstly write in AUPR scores
        f.write('AUPR score: ' + '\n\n')
        for i in range(len(method_names)):
            f.write(method_names[i] + ': ')
            f.write(str(c_AUPR[i]) + '\n')
        f.write('\n')
        #Now time to write in the missed
        for i in range(len(method_names)):
            f.write(method_names[i] + ' connection information\n\n')
            f.write('FOUND:\n')
            for j in range(int(c_found[i])):
                f.write(str(j) + ' --- ' + str(j))
                f.write('\n')
            f.write('\n')
            f.write('MISSED:\n')
            for j in range(int(c_missed[i])):
                f.write(str(j) + ' --- ' + str(j))
                f.write('\n')
            f.write('\n')
            f.write('INCORRECT:\n')
            for j in range(int(c_incorrect[i])):
                f.write(str(j) + ' --- ' + str(j))
                f.write('\n')
            f.write('\n\n')

with open('Metrics_beeline_messy.txt', 'w') as f:
        #Firstly write in AUPR scores
        f.write('AUPR score: ' + '\n\n')
        for i in range(len(method_names)):
            f.write(method_names[i] + ': ')
            f.write(str(m_AUPR[i]) + '\n')
        f.write('\n')
        #Now time to write in the missed
        for i in range(len(method_names)):
            f.write(method_names[i] + ' connection information\n\n')
            f.write('FOUND:\n')
            for j in range(int(m_found[i])):
                f.write(str(j) + ' --- ' + str(j))
                f.write('\n')
            f.write('\n')
            f.write('MISSED:\n')
            for j in range(int(m_missed[i])):
                f.write(str(j) + ' --- ' + str(j))
                f.write('\n')
            f.write('\n')
            f.write('INCORRECT:\n')
            for j in range(int(m_incorrect[i])):
                f.write(str(j) + ' --- ' + str(j))
                f.write('\n')
            f.write('\n\n')

