#This file here takes in large amounts of data, some knowledge about the most variable Transcription factors and reference network

import numpy as np
from bisect import bisect_left

#Step 1: Extract the most variable genes
#Number of Transcription factors we want
num_TFs = 100
#An empty array which holds the transcription factor
TF_interest = []
#The file we are reading in
infile = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Actual_data_tests/Research_data/inputs/scRNA-Seq/mESC/GeneOrdering.csv'
file = open(infile)

#This for loop gets the transcription factors of interest
for i, line in enumerate(file):
    if (i <= num_TFs) and (i != 0):
        #Set up an empty place for the current transcription factor
        current_TF = ''
        current = 0
        #Get it's name
        while line[current] != ',':
            current_TF = current_TF + line[current]
            current = current + 1
        #Append it
        TF_interest.append(current_TF)


#Step 2: Now get a .csv file with only the transcription factors we want to keep
#Gets the expression file
expresion_file = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Actual_data_tests/Research_data/inputs/scRNA-Seq/mESC/ExpressionData.csv'
expression = np.genfromtxt(expresion_file, delimiter= ',', dtype=str)

#Gets a name of all the transcription factors
all_TFs = expression[1:,0]
#Get the number of columns needed
num_col = len(expression[1,:])
#Holds our rows
row = []
#Find the positions of these transcription factors
for TF in TF_interest:
    #Gets the position needed
    pos = bisect_left(all_TFs,TF,0,len(all_TFs))+1
    #Adds to the rows
    row.append(expression[pos,:])

#Create the new numpy matrix
new_expression = np.array(row)

#This here saves the 
# np.savetxt("matrix_real.csv", new_expression, delimiter=",", fmt='%s')

#Step 3: Write a reference network:
reference_file = '/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Actual_data_tests/chip_x.txt'
#We now write this into a numpy array (THERE ARE 133K ROWS HERE, THIS COULD TAKE A WHILE)
reference = np.genfromtxt(reference_file, delimiter= '\t', dtype=str)

#Getting locations of the Genes
source = list(reference[1:,0])
source = [str(x) for x in source]
target = list(reference[1:,2])
target = [str(x) for x in target]

#Here I create a numpy array:
network = np.zeros(shape=(num_TFs,num_TFs))
#Convert to uppercase
TF_interest_upper = [x.upper() for x in TF_interest]

#We find all the locations
for j in range(num_TFs):
    #Firstly deal with the source
    if TF_interest_upper[j] in source:
        gene_loc_source = [i for i, x in enumerate(source) if x == TF_interest_upper[j]]
        souce_found = 1
    else:
        source_found = 0
    #We now will check every location
    if source_found == 1:
        for k in range(len(gene_loc_source)):
            #Get the target gene
            gene = target[gene_loc_source[k]]
            #See if it's a gene we needed
            if gene in TF_interest_upper:
                index = TF_interest_upper.index(gene)
                network[j][index] = 1
                network[index][j] = 1
    if TF_interest_upper[j] in target:
        gene_loc_target = [i for i, x in enumerate(target) if x == TF_interest_upper[j]]
        target_found = 1
    else:
        target_found = 0
    if target_found == 1:
        for k in range(len(gene_loc_target)):
            #Get the target gene
            gene = source[gene_loc_target[k]]
            #See if it's a gene we needed
            if gene in TF_interest_upper:
                index = TF_interest_upper.index(gene)
                network[j][index] = 1
                network[index][j] = 1







