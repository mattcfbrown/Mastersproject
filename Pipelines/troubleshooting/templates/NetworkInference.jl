
using NetworkInference

Values = ARGS


#Gets the input
genes = get_nodes(Values[1]);
network = InferredNetwork(PIDCNetworkInference(),genes);
adjacency_matrix, labels_to_ids, ids_to_labels = get_adjacency_matrix(network, parse(Float64, Values[2]))
print(adjacency_matrix)

