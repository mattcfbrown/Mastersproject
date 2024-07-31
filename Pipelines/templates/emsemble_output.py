#A way to test out the ensemble outputs, including plotting

import numpy as np
import sys
from sklearn.cluster import KMeans

# ensemble_clean = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Ensemble_testing/results/ensemble/ensemble_500.csv"
# ensemble_messy = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Ensemble_testing/results/ensemble/ensemble_500cells_messy.txt.csv"

# ensemble_clean = np.genfromtxt(ensemble_clean, delimiter=",")
# ensemble_messy = np.genfromtxt(ensemble_messy, delimiter=",")

# ensemble_clean = np.array([x for row in ensemble_clean for x in row])
# ensemble_messy = np.array([x for row in ensemble_messy for x in row])

# #Firstly let's plot what we have
# plt.figure(1)
# plt.hist(ensemble_clean ,bins = 20, histtype='bar', density = True)
# plt.title("Histogram of clean linearised points")
# plt.xlabel("Linear values")
# plt.ylabel("Density")
# plt.show()
# plt.figure(2)
# plt.hist(ensemble_messy ,bins = 20, histtype='bar', density = True)
# plt.title("Histogram of messy linearised points")
# plt.xlabel("Linear values")
# plt.ylabel("Density")
# plt.show()

# #We now should attempt to try and cluster the data we have obtained. Let's use k-means
# # Now attempt to plot the clusters with the histogram:
# kmeans_clean_final = KMeans(n_clusters=4)
# kmeans_clean_final.fit(ensemble_clean.reshape(-1,1))
# print(kmeans_clean_final.labels_)
# print(kmeans_clean_final.cluster_centers_)
# print(ensemble_clean[4])

#This here is how we can plot the data
# plt.figure(3)
# plt.scatter([i for i in range(len(ensemble_clean))], ensemble_clean, c=labels)
# plt.show()

#--------------------------------Actual code--------------------------------------

#Read in the data
input = sys.argv[1]
id = sys.argv[2]
ensemble = np.genfromtxt(input, delimiter=",")
ensemble_list = np.array([x for row in ensemble for x in row])
num_genes = int(np.size(ensemble, 1))

#Get clusters using Kmeans
kmeans = KMeans(n_clusters=4)
kmeans.fit(ensemble_list.reshape(-1,1))
max_index = np.where(kmeans.cluster_centers_ == max(kmeans.cluster_centers_ )[0])
max_index = max_index[0][0]

ensemble_matrix = np.zeros((num_genes,num_genes))
for i in range(num_genes):
    for j in range(num_genes):
        #If the label is the same as the max index, we accept it
        if (kmeans.labels_[25*i + j] == max_index):
            ensemble_matrix[i][j] = 1
        else:
            ensemble_matrix[i][j] = 0

np.savetxt("ensemble_" + str(id) + "_matrix" + ".csv", ensemble_matrix, delimiter=",")





