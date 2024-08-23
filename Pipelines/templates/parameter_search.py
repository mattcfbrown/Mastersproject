#Here I will assess the range of thresholds I can use

import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
#These two functions here serve to get me my folder
from os import listdir
from os.path import isfile, join

#I will assess two files
#25 gene file
path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/Data/outputs/grn_input_1000_clean_folder/lists'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

eb_full = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/EB/EB_output/EB_1000_clean_full.csv'
eb_ten = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/EB/EB_output/EB_1000_clean_ten.csv'
onlyfiles.append(eb_full)
onlyfiles.append(eb_ten)

#We now read in the ground truth network
original = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/EB/EB_output/truth_1000_clean.csv'
original = np.genfromtxt(original,delimiter=',')
num_genes = original.shape[0]
original[0][0] = 0

#This is the ID value
id = str(1000)

#Ensemble reading in

ensemble_full = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/Data/optimised/genci_input_1000_clean_folder/all/final_list.csv'
ensemble_ten = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/Data/optimised/genci_input_1000_clean_folder/top_10/final_list.csv'
onlyfiles.append(ensemble_full)
onlyfiles.append(ensemble_ten)

#10 gene file
# path = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/Data/outputs/grn_input_chain_clean_folder/lists'
# onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

# eb_full = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/EB/EB_output/EB_chain_clean_full.csv'
# eb_ten = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/EB/EB_output/EB_chain_clean_ten.csv'
# onlyfiles.append(eb_full)
# onlyfiles.append(eb_ten)

# #We now read in the ground truth network
# original = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/EB/EB_output/truth_chain_clean.csv'
# original = np.genfromtxt(original,delimiter=',')
# num_genes = original.shape[0]
# original[0][0] = 0

# #This is the ID value
# id = str(1000)

# #Ensemble reading in

# ensemble_full = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/Data/optimised/genci_input_chain_clean_folder/all/final_list.csv'
# ensemble_ten = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/Data/optimised/genci_input_chain_clean_folder/top_10/final_list.csv'
# onlyfiles.append(ensemble_full)
# onlyfiles.append(ensemble_ten)

#------------------------------Where all the magic happens-------------------------------
scores = [] # All the scores used
num_connections = num_genes*num_genes - num_genes

threshold = [0.05,0.1,0.2,0.3,0.4,0.5]
# threshold = np.arange(0.05, 0.5, 0.01)

for thresh in threshold:
    names = []     #Holds the names of the methods
    AUPR = []      #Holds AUPR scores
    for method in onlyfiles:
        #Firstly get the name of the method
        name = method.replace("GRN_", "")
        name = name.replace(".csv", "")
        name = name.replace("_" + id, "")
        name = name.replace("all/final_list", "ensemble_all")
        name = name.replace("top_10/final_list", "ensemble_ten")
        name = name.replace("EB_clean_ten", "EB_ten")
        name = name.replace("EB_clean_full", "EB_full")
        name = name.replace('/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/Data/optimised/genci_input_clean_folder/',"")
        name = name.replace('/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/GENCI_nextflow_pipeline/results/EB/EB_output/', "")
        names.append(name)

        # Secondly, we want to create a file 
        if name == 'EB_full':
            matrix = np.genfromtxt(eb_full, delimiter=',',
                            dtype=None)
        elif name == 'EB_ten':
            matrix = np.genfromtxt(eb_ten, delimiter=',',
                            dtype=None)
        else:
            if name == 'ensemble_all':
                file = np.genfromtxt(ensemble_full, delimiter=',',
                            dtype=None)
            elif name == 'ensemble_ten':
                file = np.genfromtxt(ensemble_ten, delimiter=',',
                            dtype=None)
            else:
                file = np.genfromtxt(path + '/' + method, delimiter=',',
                            dtype=None)
            #We only want the top 10% of connections
            matrix = np.zeros(shape=(num_genes,num_genes))
            hit = 0
            num_rows = file.shape[0]
            for i in range(len(file)):
                gene1 = int(file[i][0][1:])
                gene2 = int(file[i][1][1:])
                #Ad Hoc solution to an error I am getting
                if num_genes == 20:
                    gene1 = gene1-1
                    gene2 = gene2-1
                if matrix[gene1][gene2] == 0:
                    matrix[gene1][gene2] = 1
                    matrix[gene2][gene1] = 1
                    hit = hit + 1
                if hit == round(num_rows*thresh*0.5) or num_rows == i+1:
                    break
        
        #Now get the AUPR scores
        matrix_list = [x for row in matrix for x in row]
        original_list = [x for row in original for x in row]
        precision, recall, thresholds = metrics.precision_recall_curve(original_list,matrix_list)
        AUPR.append(metrics.auc(recall,precision))
    scores.append(AUPR)


boxplot_AUPR = {}
unique = []
for i in range(len(scores)):
    key, value = threshold[i], scores[i]
    boxplot_AUPR[key] = value
    unique.append(len(set(scores[i])))

fig, ax = plt.subplots()
ax.boxplot(boxplot_AUPR.values())
plt.xlabel("%Threshold values")
plt.ylabel("AUPR scores")
plt.title("Box plot of AUPR scores for 25 gene dataset")
ax.set_xticklabels(boxplot_AUPR.keys())
plt.show()

fig, ax = plt.subplots()
plt.plot(threshold,unique)
plt.xlabel("%Threshold values")
plt.ylabel("Number of unique values")
plt.title("Number of unique values for each threshold for 25 gene dataset")
plt.show()


