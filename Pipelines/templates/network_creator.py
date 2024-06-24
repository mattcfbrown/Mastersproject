#Here I am making a program to generate networks I am using:
#Using the networkX package

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

mydata = np.genfromtxt('/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/Data/Original_test.csv', delimiter=',')

labels = list(range(25))
labels = {labels[i]: i for i in range(25)}

def show_graph_with_labels(adjacency_matrix, mylabels):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=500, labels=mylabels, with_labels=True)
    plt.show()


show_graph_with_labels(mydata, labels)