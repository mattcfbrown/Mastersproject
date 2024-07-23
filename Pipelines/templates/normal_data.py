#This here is how we are planning on generating some normally distributed data

import pandas as pd
import numpy as np
from numpy import random as rd


#This will work in the following steps

#1: Create 100*1000 normally distributed samples
#2: Put it into tab seperated data
    #2.a: Put the numbers up the top
    #2.b: Put the gene names on the side
#Export it

data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/PIDC_testing/Data/100_ecoli1_large.txt"

rd.seed(992334)
num_cells = 2100
num_genes = 100

#All data
data = pd.read_table(data, header = None)
data = data.to_numpy()
data = data[1:,1:]
data = [x for row in data for x in row]

#We not randomise it
print(type(data))
rd.shuffle(data)

#Write it in 
f = open("shuffled.txt", "w")
f.write('0\t1000\n\n') #2.A puts numbers at the top
for i in range(0,num_genes):
    f.write('G' + str(i+1)) #2.B Put gene names on the side
    for j in range(0, num_cells):
        f.write('\t' + str(data[i+j]))
    f.write('\n')

#--------------------------------------------------------------------------------------------------
#Step 1: Getting the randomly generated data
# values = rd.normal(size= num_cells*num_genes)


#Step 2: We put it into a file
# f = open("indepenent.txt", "w")
# f.write('0\t1000\n\n') #2.A puts numbers at the top
# for i in range(0,num_genes):
#     f.write('G' + str(i+1)) #2.B Put gene names on the side
#     for j in range(0, num_cells):
#         f.write('\t' + str(values[i+j]))
#     f.write('\n')