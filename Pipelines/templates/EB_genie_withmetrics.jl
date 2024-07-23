#This file will contain the following:
#1 A calculation of Empirical Bayes with the updated prior formula
#2 Convert it into a matrix
#3 Get a list of the genes (very basic T1 to T(num_genes))
#4 Get the metric scores
#5 Get a table with the data
#6 Plot the data

#The AUPR and the gene table will be added with a python file 

#Load in all the functions
using NetworkInference
using Distributions
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

#All the values needed to perform the calculations
#This here gets the test statistics
data = Values[1]
prior_data = Values[2]
p_val = parse(Float64, Values[3])

genes = get_nodes(data);
num_cells = Int(parse(Float64, Values[4]))
if num_cells < 100
    inference = PIDCNetworkInference() #Formally MI
    distr = :Gamma
    proportion_to_keep = 0.9
else
    inference = PIDCNetworkInference() #Formally PUC
    distr = :Normal
    proportion_to_keep = 0.9
end
network_EB = InferredNetwork(inference, genes);
edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]

#Here is the prior information we will use
prior_file = DataFrame(CSV.File(prior_data, header= false))
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
        priors[(i,j)] = max(prior_file[edge1,edge2],prior_file[edge2,edge1])*2.2
    end
end
priors = [ get(priors, to_index(e.nodes), 0) for e in edge_list ]


num_bins = 5
tail = :two
w0=2.2
null_value = -Inf
multi = 10.0

function get_distr(d::Symbol)
    distrs = Dict(
        :Gamma => Gamma,
        :Normal => Normal
    )
    return distrs[d]
end

#------------------------------Part 1------------------------------
midpoints, counts, bin_width = discretize_test_statistics(test_statistics, num_bins)
null_distr = 0.5 #PIDC
# null_distr = fit_null_distribution(midpoints, counts, num_bins, bin_width, proportion_to_keep, get_distr(distr))
mixture_pdf = fit_mixture_distribution(midpoints, counts, bin_width)

#From here we will be analysing the Posterior calculation
# null_pdf(x) = pdf(null_distr, x)

function prior_fn(x) 
if x == 0*2.2
    exp(w0) / ( exp(w0) + exp(x) )
elseif 0 < x < 0.1*2.2
    exp(w0) / ( exp(w0) + 0*exp(x) )
elseif 0.1*2.2 <= x < 0.5*2.2
    exp(w0) / ( exp(w0) + multi*0.5*exp(x) )        
else
    exp(w0) / ( exp(w0) + multi*exp(x) )       
end        
end


#We hold the values here
num_test_statistics = length(test_statistics)
null_val = Array{Float64}(undef, num_test_statistics)
mix_val = Array{Float64}(undef, num_test_statistics)
posterior = Array{Float64}(undef, num_test_statistics)

#Here is where the magic happens
for i in 1:num_test_statistics
    ts = test_statistics[i]
    prior_val = prior_fn(priors[i])
    # null_val[i] = null_pdf(ts)
    null_val[i] = null_distr
    mix_val[i] = mixture_pdf(ts)

    # if mixture distr equals 0, then just return a 0 posterior
    if mix_val[i] == zero(mix_val[i])
        posterior[i] = 0.0
        continue
    end

    fdr = null_val[i] / mix_val[i]
    p1 = 1 - prior_val * fdr
    posterior[i] = p1
end

#List where the edges will be stored
eb_edges = Array{Edge}(undef, length(edge_list))
for i in 1:length(edge_list)
nodes = edge_list[i].nodes
eb_edges[i] = Edge(nodes, posterior[i])
end

#------------------------------Part 2------------------------------
eb_edges = Array{Edge}(undef, length(edge_list))
unsorted_weight = zeros(length(edge_list))
for i in 1:length(edge_list)
    nodes = edge_list[i].nodes
    eb_edges[i] = Edge(nodes, posterior[i])
    unsorted_weight[i] = eb_edges[i].weight
end

permvec = sortperm(unsorted_weight)

#Here is the final touch
matrix = zeros(Int,num_genes,num_genes)
for i in reverse(permvec)
    edge1 = eb_edges[i].nodes[1].label
    edge2 = eb_edges[i].nodes[2].label
    weight = eb_edges[i].weight
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end

#writes the file to a matrix
name = join(["EB_matrix_" ,Values[4], ".csv"])
writedlm( name,  matrix, ',')

#------------------------------Part 3------------------------------
#The variable gene_names hold all the gene names (lol)

#------------------------------Part 4------------------------------
#Here we return the Null and mixture distribution values, along with the actual data
final_data = Array{Float64}[]
for i = [null_val,mix_val,test_statistics]
    push!(final_data, i)
end

#Now put this data into a csv file with the following rows:
#Row 1 = Null values
#Row 2 = Mixture values
#Row 3 = Actual data

name_data = join(["data_" ,Values[4], ".csv"])
writedlm( name_data,  final_data, ',')

#We also wish to create a matrix where we have the data labelled, so we know where each datapoint came from 
#This loops over all the edges, and we simply save the edges in an array
gene_list = Array{String}[]
for i = 1:length(eb_edges)
    edge1 = eb_edges[i].nodes[1].label
    edge2 = eb_edges[i].nodes[2].label
    push!(gene_list, [edge1,edge2])
end

gene_name = join(["gene_list_" ,Values[4], ".txt"])
writedlm( gene_name,  gene_list, '\t')