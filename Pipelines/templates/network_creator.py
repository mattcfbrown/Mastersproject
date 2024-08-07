#Here I am making a program to generate networks I am using:
#Using the networkX package

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

mydata = np.genfromtxt('/Users/mbrown/Desktop/Research/Mastersproject/GENECI/input_data/simulated_scratch/GS/sim_eipo-modular_size-20_mixed_gs.csv', 
                       delimiter=',', dtype=None, encoding=None, skip_header=True)

list_of_lists = []
for i in range(mydata.size):
    temp = []
    for j in range(len(mydata[i])-1):
        temp.append(mydata[i][j+1])
    list_of_lists.append(temp)

matrix = np.array(list_of_lists)

# gene_dict = {genes[i]: genes[i] for i in range(num_genes)}
# matrix = np.zeros((num_genes,num_genes))
num_genes = matrix[0].size
labels = list(range(num_genes))
labels = {labels[i]: i for i in range(num_genes)}

def show_graph_with_labels(adjacency_matrix, mylabels):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.DiGraph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=500, labels=mylabels, with_labels=True)
    plt.show()

rows, cols = np.where(matrix == 1)
edges = zip(rows.tolist(), cols.tolist())
gr = nx.DiGraph()
gr.add_edges_from(edges)
print(nx.find_cycle(gr))
print(nx.simple_cycles(gr))

show_graph_with_labels(matrix, labels)