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

Values = ARGS
data = Values[1]
w0 = parse(Float64, Values[2])
full_data = Values[3]
genie_data = Values[4]
zero_data = Values[5]
p_val = parse(Float64, Values[6])
id = Values[7]

dist = :Normal
inference = PUCNetworkInference()

#Gets the network
genes = get_nodes(data);
network_EB = InferredNetwork(inference, genes);

#Prior files
full_file = DataFrame(CSV.File(full_data, header= false))
full_file = Matrix(full_file)
genie_file = DataFrame(CSV.File(genie_data, header= false))
genie_file = Matrix(genie_file)
zero_file = DataFrame(CSV.File(zero_data, header= false))
zero_file = Matrix(zero_file)

#Now generates a dictionary of the prior information
num_genes = length(genes)
#Generate the names of the gene
gene_names = []
for i in 0:(num_genes-1)
    cur_gene = string("T", i)
    push!(gene_names,cur_gene)
end
#Builds the dictionary
full_priors = Dict()
genie_priors = Dict()
zero_priors = Dict()
for i in gene_names
    for j in gene_names
        edge1 = parse(Int,i[2:end]) + 1
        edge2 = parse(Int,j[2:end]) + 1
        full_priors[(i,j)] = full_file[edge1,edge2]*2.2
        genie_priors[(i,j)] = genie_file[edge1,edge2]*2.2
        zero_priors[(i,j)] = zero_file[edge1,edge2]*2.2
    end
end

edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]
full_list  = [ get(full_priors, to_index(e.nodes), 0) for e in edge_list ]
genie_list = [ get(genie_priors, to_index(e.nodes), 0) for e in edge_list ]
zero_list  = [ get(zero_priors, to_index(e.nodes), 0) for e in edge_list ]

#Now run Empirical Bayes
num_bins = 5
distr = dist
proportion_to_keep = 0.8
if id == "100cells_messy"
    proportion_to_keep = 0.6
end
if id == "100cells"
    proportion_to_keep = 0.6
end
tail = :two
full_posteriors = empirical_bayes(test_statistics, full_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)
genie_posteriors = empirical_bayes(test_statistics, genie_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)
zero_posteriors = empirical_bayes(test_statistics, zero_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)

#-------Full-------
#List to hold edges
edges = Array{Edge}(undef, length(edge_list))
#We also want to get weighted values
weighted = zeros(length(edge_list))
#Gets the edges
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes

    edges[i] = Edge(nodes, full_posteriors[i])
    weighted[i] = edges[i].weight
end

permvec = sortperm(weighted)
#We now want to create a matrix which holds all the values
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

#writes the file to a matrix
name = join(["Eb_fullmatrix_", id, "_", Values[2], ".csv"])
writedlm( name,  matrix, ',')


#-----GENIE--------
#List to hold edges
edges = Array{Edge}(undef, length(edge_list))
#We also want to get weighted values
weighted = zeros(length(edge_list))
#Gets the edges
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes

    edges[i] = Edge(nodes, genie_posteriors[i])
    weighted[i] = edges[i].weight
end

permvec = sortperm(weighted)
#We now want to create a matrix which holds all the values
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
#writes the file to a matrix
name = join(["Eb_geniematrix_", id, "_", Values[2], ".csv"])
writedlm( name,  matrix, ',')

#------Zero--------
#List to hold edges
edges = Array{Edge}(undef, length(edge_list))
#We also want to get weighted values
weighted = zeros(length(edge_list))
#Gets the edges
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes

    edges[i] = Edge(nodes, zero_posteriors[i])
    weighted[i] = edges[i].weight
end

permvec = sortperm(weighted)
#We now want to create a matrix which holds all the values
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
#writes the file to a matrix
name = join(["Eb_zeromatrix_", id, "_", Values[2] ,".csv"])
writedlm( name,  matrix, ',')
