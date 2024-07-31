#This will attempt to create fancy stutter box plots for each scenario to determine
#Normalisation effectivness

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys


names = ['Genie', 'tigress', 'pearson', 'CLR'] #The names of all the ensemble cklassifiers
#Now get the original
original = np.genfromtxt(sys.argv[1], delimiter=",")
original[0][0] = 0
num_genes = int(np.size(original, 1))
original = [x for row in original for x in row]
#This is the id (the number of cells)
id = sys.argv[2]

#We perform the following loop for each plots
for i in range(7,len(sys.argv)):
    scores = np.genfromtxt(sys.argv[i],delimiter=",")
    scores = [x for row in scores for x in row]         #Normalisation scores
    guess = np.genfromtxt(sys.argv[i-4],delimiter=",") 
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
    
    #Now it is time to plot it
    plot_name = 'normalisation_' + str(id) + '_' + names[i-7] + '.pdf'
    # plt.boxplot(A_correctly + A_incorrect + R_correctly + R_incorrect, whis=100, widths = 0)
    # plt.scatter(np.random.normal(1, 0.04, size=len(A_correctly)), A_correctly, color = 'b')
    # plt.scatter(np.random.normal(1, 0.04, size=len(A_incorrect)), A_incorrect, color = 'r')
    # plt.scatter(np.random.normal(1, 0.04, size=len(R_correctly)), R_correctly, color = 'g')
    # plt.scatter(np.random.normal(1, 0.04, size=len(R_incorrect)), R_incorrect, color = 'y')
    # A_C = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
    #                       markersize=10, label='Accepted correctly')
    # A_I = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
    #                       markersize=10, label='Accepted Incorrectly')
    # R_C = mlines.Line2D([], [], color='green', marker='o', linestyle='None',
    #                       markersize=10, label='Rejected correctly')
    # R_I = mlines.Line2D([], [], color='yellow', marker='o', linestyle='None',
    #                       markersize=10, label='Rejected Incorrectly')
    # plt.legend(handles=[A_C,A_I,R_C,R_I])
    # plt.title("Plot of linearlised points")
    # plt.xlabel(" ")
    # plt.ylabel("Linearlised value")
    plt.figure(i-6)
    plt.hist((A_correctly,A_incorrect,R_correctly,R_incorrect), bins = 20, histtype='bar', stacked=True, density = True)
    plt.legend(['Accepted correctly','Accepted incorrectly', 'Rejected correctly', 'Rejected incorrectly'])
    plt.title("Histogram of linearised points")
    plt.xlabel("Linear values")
    plt.ylabel("Density")
    plt.savefig(plot_name, format = "pdf") 
    plt.show()   