#Here I aim to plot and see how valuable the GENIE normalisation techniques are

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import re
import math

#File:
infile = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/PIDC_testing/results/geniemetrics/genie_outputs/genie1000.txt"
truth = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/PIDC_testing/Data/Original_test.csv"
truth = np.genfromtxt(truth,delimiter=',')
num_genes = 25

#This produces a list containing all the values
data = []
with open(infile, 'r') as fp:
    lines = fp.readlines()

for line in lines[1:]:
    #Gets the location of the genes
    location =[loc.start() for loc in re.finditer('Gene', line)]
    #This gives us the gene number
    location = [x+4 for x in location]
    location.append(0)
    vals = []
    for loc in location:
        current = loc
        gene_ID = ''
        while line[current] != ' ' and current < len(line)-1:
            gene_ID = gene_ID + line[current]
            current = current + 1
        location[2] = location[1] + len(gene_ID) + 1   
        vals.append(gene_ID)
    data.append(vals)


#This is the matrix stuff
matrix = np.zeros((num_genes,num_genes))

#Min and max of the data
max = float(data[0][2])
min = float(data[math.perm(num_genes,2)-1][2])

#--------------------------Plotting-----------------------------------
genie_points = []    #These are all the GENIE values
actual_data = []     #This is the data inputs which actually are found
for value in data:
    if truth[int(value[0])-1][int(value[1])-1] == 1:
        actual_data.append(float(value[2]))
    else:      
        genie_points.append(float(value[2])) 

#Firstly a histogram of the Genie values
plt.figure(1)
plt.hist((genie_points,actual_data), bins = 20, histtype='bar', stacked=True, density = True)
plt.legend(['No connection','Connection'])
plt.title("Histogram of Genie points")
plt.xlabel("GENIE values")
plt.ylabel("Density")
plt.axvline(0.10,color = 'r')

#Secondly a plot with linear transformed data
#This is for the genie points not found
genie_linear = []
for i in range(len(genie_points)):
    numerator = genie_points[i]-min
    denominator = max-min
    genie_linear.append(numerator/denominator)

#Now all true points
actual_linear = []
for i in range(len(actual_data)):
    numerator = actual_data[i]-min
    denominator = max-min
    actual_linear.append(numerator/denominator)

plt.figure(2)
plt.boxplot(genie_linear + actual_linear, whis=100, widths = 0)
plt.scatter(np.random.normal(1, 0.04, size=len(genie_linear)), genie_linear, color = 'b')
plt.scatter(np.random.normal(1, 0.04, size=len(actual_linear)), actual_linear, color = 'r')
no_connection = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
                          markersize=10, label='No connections')
Connection = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
                          markersize=10, label='Connections')
plt.legend(handles=[no_connection,Connection])
plt.title("Plot of Genie linearlised points")
plt.xlabel(" ")
plt.ylabel("Linearlised value")

#Now we will normalise it:
mu = sum(genie_points + actual_data)/len(genie_points + actual_data)
sd = np.std(genie_points + actual_data)
def normalised(x,mu,sd):
    return (x-mu)/sd
genie_normalised = [normalised(x,sd,mu) for x in genie_points]
actual_normalised = [normalised(x,sd,mu) for x in actual_data]

plt.figure(3)
plt.boxplot(genie_normalised + actual_normalised, whis=100, widths = 0)
plt.scatter(np.random.normal(1, 0.04, size=len(genie_normalised)), genie_normalised, color = 'b')
plt.scatter(np.random.normal(1, 0.04, size=len(actual_normalised)), actual_normalised, color = 'r')
no_connection = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
                          markersize=10, label='No connections')
Connection = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
                          markersize=10, label='Connections')
plt.legend(handles=[no_connection,Connection])
plt.title("Plot of Genie normalised points")
plt.xlabel(" ")
plt.ylabel("Normalised value")


plt.show()
input()
