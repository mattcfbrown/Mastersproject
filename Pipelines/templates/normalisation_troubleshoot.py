#This file will be deleted later, but it intends to test how good the normalisation visulisation plotting is

import numpy as np
import matplotlib.pyplot as plt

#We have three madeup plots:

scores = np.array([[1,10,1],[10,10,1],[10,1,10]]) #Scores
scores = [x for row in scores for x in row]
original = np.array([[0,1,1],[1,0,0],[0,1,0]])    #Original
original = [x for row in original for x in row]
guess = np.array([[0,1,0],[1,0,1],[0,1,0]])       #Guess
guess = [x for row in guess for x in row]
index = range(1,10)


A_correctly = [] #Accepted correctly
A_incorrect = [] #Accepted incorrectly
R_correctly = [] #Rejected correctly
R_incorrect = [] #Rejected incorrectly
for i in range(0,9):
    if original[i] == 1:
        #Correctly
        if guess[i] == 1:
            A_correctly.append(scores[i])
        #Incorrectly
        else:
            R_incorrect.append(scores[i])
    #Rejected
    else:
        #Correctly
        if guess[i] == 0:
            R_correctly.append(scores[i])
        #Incorrectly
        else:
            A_incorrect.append(scores[i])

plt.figure(1)
plt.hist((A_correctly,A_incorrect,R_correctly,R_incorrect), bins = 20, histtype='bar', stacked=True, density = True)
plt.legend(['Accepted correctly','Accepted incorrectly', 'Rejected correctly', 'Rejected incorrectly'])
plt.title("Histogram of linearised points")
plt.xlabel("Linear values")
plt.ylabel("Density")
plt.show()   