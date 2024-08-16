#Here I will run Empirical Bayes using GENECI prior file

using NetworkInference
using EmpiricalBayes
using DelimitedFiles
using DataFrames, CSV

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

#Input and important values
Values = ARGS
data = Values[1] #This is the data
#We now read in the input
genes = get_nodes(data);
inference = PUCNetworkInference()                #We wish to use PUC for this
network_EB = InferredNetwork(inference, genes);
#Get the number of genes
num_genes = length(genes)


#Now we get the priors
full_prior = Values[2]
full_prior = CSV.read(full_prior, DataFrame, header=false)
#We convert this into a matrix
full_matrix = Dict()
for i in 1:nrow(full_prior)
    edge1 = parse(Int,full_prior[i,1][2:end]) + 1
    edge1_name = join(["T", string(edge1)])
    edge2 = parse(Int,full_prior[i,2][2:end]) + 1
    edge2_name = join(["T", string(edge2)])
    full_matrix[(edge1_name,edge2_name)] = full_prior[i,3]
end

#For the 10 priors
ten_prior = Values[3]
ten_prior = CSV.read(ten_prior, DataFrame, header=false)
#We convert this into a matrix
ten_matrix = Dict()
for i in 1:nrow(ten_prior)
    edge1 = parse(Int,ten_prior[i,1][2:end]) + 1
    edge1_name = join(["T", string(edge1)])
    edge2 = parse(Int,ten_prior[i,2][2:end]) + 1
    edge2_name = join(["T", string(edge2)])
    ten_matrix[(edge1_name,edge2_name)] = ten_prior[i,3]
end

#Edge list
edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]

#Get the prior lists
full_list = [ get(full_matrix, to_index(e.nodes), 0) for e in edge_list ]
ten_list = [ get(ten_matrix, to_index(e.nodes), 0) for e in edge_list ]

#Some variables
num_bins = 5
distr = :Normal
proportion_to_keep = 0.9
tail = :two
w0 = 2.2

posteriors_full = empirical_bayes(test_statistics, full_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)
posteriors_ten = empirical_bayes(test_statistics, ten_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)


#List to hold edges
edges = Array{Edge}(undef, length(edge_list))
#We also want to get weighted values
weighted = zeros(length(edge_list))
#Gets the edges
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes

    edges[i] = Edge(nodes, posteriors_full[i])
    weighted[i] = edges[i].weight
end

p_val = parse(Float64,Values[4])

permvec = sortperm(weighted)
matrix = zeros(Int,num_genes,num_genes)
for i in reverse(permvec)
    edge1 = edges[i].nodes[1].label
    edge2 = edges[i].nodes[2].label
    weight = edges[i].weight
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end
num_cells = Values[5]
name = join(["EB_", num_cells, "_full.csv"])
writedlm( name,  matrix, ',')


#----------I should just write a function but I am lazy------------
#List to hold edges
edges = Array{Edge}(undef, length(edge_list))
#We also want to get weighted values
weighted = zeros(length(edge_list))
#Gets the edges
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes

    edges[i] = Edge(nodes, posteriors_ten[i])
    weighted[i] = edges[i].weight
end

p_val = parse(Float64,Values[4])

permvec = sortperm(weighted)
matrix = zeros(Int,num_genes,num_genes)
for i in reverse(permvec)
    edge1 = edges[i].nodes[1].label
    edge2 = edges[i].nodes[2].label
    weight = edges[i].weight
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end
num_cells = Values[5]
name = join(["EB_", num_cells, "_ten.csv"])
writedlm( name,  matrix, ',')