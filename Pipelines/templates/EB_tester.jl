using NetworkInference
using DelimitedFiles
using DataFrames, CSV


data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/results/full_priors/formatted_data.txt"
p_val = 0.9
prior_data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/Data/Original_test.csv"
type = "Test"
to_keep = 0.9
dist = :Normal

genes = get_nodes(data);
network_EB = InferredNetwork(PUCNetworkInference(), genes);

#Gets the prior information
prior_file = DataFrame(CSV.File(prior_data, header= false))
prior_file = Matrix(prior_file)


#Now generates a dictionary of the prior information
num_genes = length(genes)
#Generate the names of the gene
gene_names = []
for i in 0:(num_genes-1)
    cur_gene = string("T", i)
    push!(gene_names,cur_gene)
end
#Builds the dictionary
priors = Dict()
for i in gene_names
    for j in gene_names
        edge1 = parse(Int,i[2:end]) + 1
        edge2 = parse(Int,j[2:end]) + 1
        priors[(i,j)] = prior_file[edge1,edge2]
    end
end


edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]
prior_list = [ get(priors, to_index(e.nodes), 0) for e in edge_list ]

#List where the edges will be stored
eb_edges = Array{Edge}(undef, length(edge_list))

#Now run Empirical Bayes
num_bins = 5
distr = dist
proportion_to_keep = to_keep
tail = :two
w0 = 2.2
posteriors = empirical_bayes(test_statistics, prior_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)

for i in 1:length(edge_list)
    nodes = edge_list[i].nodes
    eb_edges[i] = Edge(nodes, posteriors[i])
end


#Here I am going to print out the weights
print(length(eb_edges))
for i =1:length(eb_edges)
    weight = eb_edges[i].weight
    print(weight)
    print("\n")
end