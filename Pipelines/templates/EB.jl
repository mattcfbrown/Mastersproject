using NetworkInference
using EmpiricalBayes
using DelimitedFiles

#This function is used in the Empirical Bayes framework
function to_index(label1::AbstractString, label2::AbstractString)
    if label1 > label2
        return (label1, label2)
    else
        return (label2, label1)
    end
end

function to_index(node1::Node, node2::Node)
    return to_index(node1.label, node2.label)
end

function to_index(nodes::Array{Node, 1})
    n1, n2 = nodes
    return to_index(n1.label, n2.label)
end


#Values = ARGS

#Gets the network
genes = get_nodes("/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/results/Empirical_Bayes/formatted_data.txt");
network_EB = InferredNetwork(MINetworkInference(), genes);

#Now generates a dictionary of the prior information
num_genes = length(genes)
#Generate the names of the gene
gene_names = []
for i in 0:(num_genes-1)
    name = string("T", i)
    push!(gene_names,name)
end
#Builds the dictionary
priors = Dict()
for i in gene_names
    for j in gene_names
        priors[(i,j)] = 0
    end
end


edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]
prior_list = [ get(priors, to_index(e.nodes), 0) for e in edge_list ]

#List where the edges will be stored
eb_edges = Array{Edge}(undef, length(edge_list))

#Now run Empirical Bayes
num_bins = 5
distr = :Gamma
proportion_to_keep = 1.0
tail = :two
w0 = 2.2
posteriors = empirical_bayes(test_statistics, prior_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)

for i in 1:length(edge_list)
    nodes = edge_list[i].nodes
    eb_edges[i] = Edge(nodes, posteriors[i])
end

# Remove infinite values
eb_edges = filter(x->isfinite(x.weight), eb_edges)
#sorts by highest to lowest p value
sort!(eb_edges, rev = true, by = x->x.weight)


#We now want to create a matrix which holds all the values
matrix = zeros(Int,num_genes,num_genes)
p_value = 0.9
for i = 1:length(eb_edges)
    edge1 = eb_edges[i].nodes[1].label
    edge2 = eb_edges[i].nodes[2].label
    weight = eb_edges[i].weight
    if weight < p_value
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end

#writes the file to a matrix
writedlm( "Eb_test.csv",  matrix, ',')