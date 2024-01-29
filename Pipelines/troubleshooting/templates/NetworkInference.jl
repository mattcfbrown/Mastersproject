
using NetworkInference

Values = ARGS


#Gets the input
genes = get_nodes(Values[1]);
threshold = Values[2]
network = InferredNetwork(PIDCNetworkInference(),genes);

#open(Values[3], "w") do file
#    write(file,network)
#end

adjacency_matrix, labels_to_ids, ids_to_labels = get_adjacency_matrix(network, parse(Float64, threshold))
display(adjacency_matrix)

#number_of_edges = false ?
#    findfirst(x -> x.weight < threshold, network.edges) - 1 :
#    Int(round(length(network.edges)*parse(Float64, threshold)))


#    Let's extract these edges:
#print("\n")
#for edge in network.edges[1 : number_of_edges]
#    node1 = labels_to_ids[edge.nodes[1].label]
#    node2 = labels_to_ids[edge.nodes[2].label]
#    print(ids_to_labels[node1], " ", ids_to_labels[node2], "\n")
#end

