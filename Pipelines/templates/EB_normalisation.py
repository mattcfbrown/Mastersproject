import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys


#Now get the original
original = np.genfromtxt(sys.argv[1], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))
original = [x for row in original for x in row]
#This is the id (the number of cells)
id = sys.argv[2]


scores = np.genfromtxt(sys.argv[3],delimiter=",")
scores = [x for row in scores for x in row]         #Normalisation scores
guess = np.genfromtxt(sys.argv[4],delimiter=",") 
guess = [x for row in guess for x in row]           #The guess
A_correctly = [] #Accepted correctly
A_incorrect = [] #Accepted incorrectly
R_correctly = [] #Rejected correctly
R_incorrect = [] #Rejected incorrectly
for j in range(len(scores)):
    #Accepted
    if original[j] == 1:
        #Correctly
        if guess[j] == 1:
            A_correctly.append(scores[j])
        #Incorrectly
        else:
            R_incorrect.append(scores[j])
    #Rejected
    else:
        #Correctly
        if guess[j] == 0:
            R_correctly.append(scores[j])
        #Incorrectly
        else:
            A_incorrect.append(scores[j])

plot_name = 'normalisation_' + str(id) + '_EB.pdf'
plt.figure(1)
plt.hist((A_correctly,A_incorrect,R_correctly,R_incorrect), bins = 20, histtype='bar', stacked=True, density = True)
plt.legend(['Accepted correctly','Accepted incorrectly', 'Rejected correctly', 'Rejected incorrectly'])
plt.title("Histogram of linearised points")
plt.xlabel("Linear values")
plt.ylabel("Density")
plt.savefig(plot_name, format = "pdf") 
plt.show()   
