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

# data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/Data/25genes_2000cells.txt"
# p_val = 0.9
# prior_data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/Data/Original_test.csv"
# type = "Test"
# to_keep = 0.9
# dist = :Normal


data = Values[1]
p_val = parse(Float64, Values[2])
prior_data = Values[3]
type = Values[4]
to_keep = parse(Float64, Values[5])
if Values[6] == "gam"
    dist = :Normal
elseif Values[6] == "normal"
    dist = :Normal
end
if Values[7] == "PUC"
    inference = PUCNetworkInference()
elseif Values[7] == "MI"
    inference = PUCNetworkInference()
end
# w0 = parse(Float64, Values[8])
w0 = 2.2


#Gets the network
genes = get_nodes(data);
network_EB = InferredNetwork(inference, genes);

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
        priors[(i,j)] = prior_file[edge1,edge2]*2.2
    end
end


edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]
prior_list = [ get(priors, to_index(e.nodes), 0) for e in edge_list ]



#Now run Empirical Bayes
num_bins = 5
distr = dist
proportion_to_keep = to_keep
tail = :two
posteriors = empirical_bayes(test_statistics, prior_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)
print(posteriors)

#List to hold edges
edges = Array{Edge}(undef, length(edge_list))
#We also want to get weighted values
weighted = zeros(length(edge_list))
#Gets the edges
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes

    edges[i] = Edge(nodes, posteriors[i])
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
name = join(["Eb_matrix_", type, ".csv"])
writedlm( name,  matrix, ',')