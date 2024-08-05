#Here I will run Empirical Bayes the required amount of times needed

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
data = Values[1] #This is the data
#We now read in the input
genes = get_nodes(data);
inference = PUCNetworkInference()                #We wish to use PUC for this
network_EB = InferredNetwork(inference, genes);
#Get the number of genes
num_genes = length(genes)
#All the priors
pearson_prior = Values[2]
spear_prior   = Values[3]
MI_prior      = Values[4]
#Converts all the prior information
pearson_file = DataFrame(CSV.File(pearson_prior, header= false))
pearson_file = Matrix(pearson_file)
spear_file   = DataFrame(CSV.File(spear_prior, header= false))
spear_file   = Matrix(spear_file)
MI_file      = DataFrame(CSV.File(MI_prior, header= false))
MI_file      = Matrix(MI_file)
#Get the zero prior
zero_file    = zeros(num_genes,num_genes)

#get the name of the genes
gene_names = []
for i in 0:(num_genes-1)
    cur_gene = string("T", i)
    push!(gene_names,cur_gene)
end

#Build the priors
pearson_priors = Dict()
spear_priors   = Dict()
MI_priors      = Dict()
zero_priors    = Dict()
for i in gene_names
    for j in gene_names
        edge1 = parse(Int,i[2:end]) + 1
        edge2 = parse(Int,j[2:end]) + 1
        pearson_priors[(i,j)] = pearson_file[edge1,edge2]*2.2
        spear_priors[(i,j)]   = spear_file[edge1,edge2]*2.2
        MI_priors[(i,j)]      = MI_file[edge1,edge2]*2.2
        zero_priors[(i,j)]    = zero_file[edge1,edge2]*2.2
    end
end

#We now get the test statistics needed
edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]

#Get the prior lists
pearson_list = [ get(pearson_priors, to_index(e.nodes), 0) for e in edge_list ]
spear_list   = [ get(spear_priors, to_index(e.nodes), 0) for e in edge_list ]
MI_list      = [ get(MI_priors, to_index(e.nodes), 0) for e in edge_list ]
zero_list    = [ get(zero_priors, to_index(e.nodes), 0) for e in edge_list ]

#Some variables
num_bins = 5
distr = :Normal
proportion_to_keep = 0.9
tail = :two
w0 = 2.2

#We now run EB
pearson_posteriors = empirical_bayes(test_statistics, pearson_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)
spear_posteriors   = empirical_bayes(test_statistics, spear_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)
MI_posteriors      = empirical_bayes(test_statistics, MI_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)
zero_posteriors    = empirical_bayes(test_statistics, zero_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)

#List to hold edges
pearson_edges = Array{Edge}(undef, length(edge_list))
spear_edges   = Array{Edge}(undef, length(edge_list))
MI_edges      = Array{Edge}(undef, length(edge_list))
zero_edges    = Array{Edge}(undef, length(edge_list))
#We also want to get weighted values
pearson_weight = zeros(length(edge_list))
spear_weight   = zeros(length(edge_list))
MI_weight      = zeros(length(edge_list))
zero_weight    = zeros(length(edge_list))
#Gets the edges
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes

    pearson_edges[i] = Edge(nodes, pearson_posteriors[i])
    pearson_weight[i] = pearson_edges[i].weight

    spear_edges[i] = Edge(nodes, spear_posteriors[i])
    spear_weight[i] = spear_edges[i].weight

    MI_edges[i] = Edge(nodes, MI_posteriors[i])
    MI_weight[i] = MI_edges[i].weight

    zero_edges[i] = Edge(nodes, zero_posteriors[i])
    zero_weight[i] = zero_edges[i].weight
end

p_val = parse(Float64,Values[5])

permvec = sortperm(pearson_weight)
matrix = zeros(Int,num_genes,num_genes)
for i in reverse(permvec)
    edge1 = pearson_edges[i].nodes[1].label
    edge2 = pearson_edges[i].nodes[2].label
    weight = pearson_edges[i].weight
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end
num_cells = Values[6]
name = join(["Eb_pearson_", num_cells, ".csv"])
writedlm( name,  matrix, ',')

permvec = sortperm(spear_weight)
matrix = zeros(Int,num_genes,num_genes)
for i in reverse(permvec)
    edge1 = spear_edges[i].nodes[1].label
    edge2 = spear_edges[i].nodes[2].label
    weight = spear_edges[i].weight
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end
name = join(["Eb_spear_", num_cells, ".csv"])
writedlm( name,  matrix, ',')


permvec = sortperm(MI_weight)
matrix = zeros(Int,num_genes,num_genes)
for i in reverse(permvec)
    edge1 = MI_edges[i].nodes[1].label
    edge2 = MI_edges[i].nodes[2].label
    weight = MI_edges[i].weight
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end
name = join(["Eb_MI_", num_cells, ".csv"])
writedlm( name,  matrix, ',')


permvec = sortperm(zero_weight)
matrix = zeros(Int,num_genes,num_genes)
for i in reverse(permvec)
    edge1 = zero_edges[i].nodes[1].label
    edge2 = zero_edges[i].nodes[2].label
    weight = zero_edges[i].weight
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end
name = join(["Eb_zero_", num_cells, ".csv"])
writedlm( name,  matrix, ',')



