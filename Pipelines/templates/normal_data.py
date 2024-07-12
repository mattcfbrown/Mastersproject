#This here is how we are planning on generating some normally distributed data

from numpy import random as rd

#This will work in the following steps

#1: Create 100*1000 normally distributed samples
#2: Put it into tab seperated data
    #2.a: Put the numbers up the top
    #2.b: Put the gene names on the side
#Export it

rd.seed(992334)
num_cells = 1000
num_genes = 100

#Step 1: Getting the randomly generated data
values = rd.normal(size= num_cells*num_genes)

#Step 2: We put it into a file
f = open("indepenent.txt", "w")
f.write('0\t1000\n\n') #2.A puts numbers at the top
for i in range(0,num_genes):
    f.write('G' + str(i+1)) #2.B Put gene names on the side
    for j in range(0, num_cells):
        f.write('\t' + str(values[i+j]))
    f.write('\n')