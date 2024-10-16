#Because I am an idiot, here is another PIDC function

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

#Inputted values:
data = Values[1]
inference = PIDCNetworkInference() #This will remain constant
num_bins = 5
w0 = 2.2
priors = Values[2]
p_val = parse(Float64,Values[3])
id = Values[4]
type = Values[5]

genes = get_nodes(data);
num_genes = length(genes)
network_EB = InferredNetwork(inference, genes);

edge_list = network_EB.edges;
# #This here will spit out the test statistics needed
test_statistics = [e.weight for e in edge_list]

gene_names = []
for i in 0:(num_genes-1)
    cur_gene = string("T", i)
    push!(gene_names,cur_gene)
end

prior = DataFrame(CSV.File(priors, header= false))
prior = Matrix(prior)
prior_matrix = Dict()
for i in gene_names
    for j in gene_names
        edge1 = parse(Int,i[2:end]) + 1
        edge2 = parse(Int,j[2:end]) + 1
        prior_matrix[(i,j)] = prior[edge1,edge2]*2.2
    end
end

prior_list  = [ get(prior_matrix, to_index(e.nodes), 0) for e in edge_list ]

#----------------------------------Where the magic happens--------------------------------------
midpoints, counts, bin_width = discretize_test_statistics(test_statistics, num_bins)
null_distr = 0.5 #Uniform distribution
mixture_pdf = fit_mixture_distribution(midpoints, counts, bin_width) #Mixture is important

prior_fn(x) = exp(w0) / ( exp(w0) + exp(x) )

num_test_statistics = length(test_statistics)
posteriors = Array{Float64}(undef, num_test_statistics)


for i in 1:num_test_statistics
    ts = test_statistics[i]
    prior_val = prior_fn(prior_list[i])
    null_val = null_distr
    mix_val = mixture_pdf(ts)

    # if mixture distr equals 0, then just return a 0 posterior
    if mix_val == zero(mix_val)
        posteriors[i] = 0.0
        continue
    end

    fdr = null_val / mix_val
    posteriors[i] = 1 - prior_val * fdr
end

#----------------------------------Formatting time-------------------------------------------
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
name = join(["Eb_", id, "_", type, ".csv"])
writedlm( name,  matrix, ',')