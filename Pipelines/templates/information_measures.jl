#Here I am going to come out with three different methods

using NetworkInference
using EmpiricalBayes
using DelimitedFiles
using DataFrames, CSV

#Inputs
Values = ARGS
data = Values[1]
id = Values[2]
genes = get_nodes(data);
num_genes = length(genes)
threshold = 0.15

#Mutual information
# network_MI = InferredNetwork(MINetworkInference(),genes); #Get the network
# weights_MI = [e.weight for e in network_MI.edges]  #Gets the values from the network
# max_mi = maximum(weights_MI) #Gets the maximum
# min_mi = minimum(weights_MI) #Gets the minimum
# MI_scores = zeros(Float64,num_genes,num_genes) #This will have the weights
# for edge in network_MI.edges
#     edge1 = parse(Int,edge.nodes[1].label[2:end]) + 1
#     edge2 = parse(Int,edge.nodes[2].label[2:end]) + 1
#     weight = edge.weight
#     MI_scores[edge1,edge2] = (weight-min_mi)/(max_mi-min_mi)
#     MI_scores[edge2,edge1] = (weight-min_mi)/(max_mi-min_mi)
# end
# name_score = join(["MI_scores_", id, ".csv"])
# writedlm(name_score,  MI_scores, ',')

# adjacency_matrix_MI, labels_to_ids, ids_to_labels = get_adjacency_matrix(network_MI, threshold);
# name_matrix = join(["MI_matrix_", id, ".csv"])
# writedlm(name_matrix,  Int.(adjacency_matrix_MI), ',')


# #PUC
# network_PUC = InferredNetwork(PUCNetworkInference(),genes); #Get the network
# weights_PUC = [e.weight for e in network_PUC.edges]  #Gets the values from the network
# max_puc = maximum(weights_PUC) #Gets the maximum
# min_puc = minimum(weights_PUC) #Gets the minimum
# PUC_scores = zeros(Float64,num_genes,num_genes) #This will have the weights
# for edge in network_PUC.edges
#     edge1 = parse(Int,edge.nodes[1].label[2:end]) + 1
#     edge2 = parse(Int,edge.nodes[2].label[2:end]) + 1
#     weight = edge.weight
#     PUC_scores[edge1,edge2] = (weight-min_puc)/(max_puc-min_puc)
#     PUC_scores[edge2,edge1] = (weight-min_puc)/(max_puc-min_puc)
# end
# name_score = join(["PUC_scores_", id, ".csv"])
# writedlm(name_score,  PUC_scores, ',')

# adjacency_matrix_PUC, labels_to_ids, ids_to_labels = get_adjacency_matrix(network_PUC, threshold);
# name_matrix = join(["PUC_matrix_", id, ".csv"])
# writedlm(name_matrix,  Int.(adjacency_matrix_PUC), ',')

#CLR
network_CLR = InferredNetwork(CLRNetworkInference(),genes); #Get the network
weights_CLR = [e.weight for e in network_CLR.edges]  #Gets the values from the network
max_clr = maximum(weights_CLR) #Gets the maximum
min_clr = minimum(weights_CLR) #Gets the minimum
CLR_scores = zeros(Float64,num_genes,num_genes) #This will have the weights
for edge in network_CLR.edges
    edge1 = parse(Int,edge.nodes[1].label[2:end]) + 1
    edge2 = parse(Int,edge.nodes[2].label[2:end]) + 1
    weight = edge.weight
    CLR_scores[edge1,edge2] = (weight-min_clr)/(max_clr-min_clr)
    CLR_scores[edge2,edge1] = (weight-min_clr)/(max_clr-min_clr)
end
name_score = join(["CLR_scores_", id, ".csv"])
writedlm(name_score,  CLR_scores, ',')

adjacency_matrix_CLR, labels_to_ids, ids_to_labels = get_adjacency_matrix(network_CLR, threshold);
name_matrix = join(["CLR_matrix_", id, ".csv"])
writedlm(name_matrix,  Int.(adjacency_matrix_CLR), ',')

