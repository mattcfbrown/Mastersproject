#Here we will run PIDC method for empirical Bayes

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
num_bins = 5 #Idk how many to do yet, this should do for now (TEST THIS IF THE METHOD WORKS)
prior_data = Values[2]
w0 = 2.2
multi = 10
p_val = parse(Float64,Values[3])

genes = get_nodes(data);
network_EB = InferredNetwork(inference, genes);
edge_list = network_EB.edges;
# #This here will spit out the test statistics needed
test_statistics = [e.weight for e in edge_list]

#------------------------------------------Priors----------------------------------------------

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
prior_list = [ get(priors, to_index(e.nodes), 0) for e in edge_list ]


#----------------------------------Where the magic happens--------------------------------------
midpoints, counts, bin_width = discretize_test_statistics(test_statistics, num_bins)
null_distr = 0.5 #Uniform distribution
mixture_pdf = fit_mixture_distribution(midpoints, counts, bin_width) #Mixture is important


#--------------------------------Method to test now-------------------------------------------

# Prior function
function prior_fn(x) 
    if x == 0
        exp(w0) / ( exp(w0) + exp(x) )
    elseif 0 < x < 0.1
        exp(w0) / ( exp(w0) + 0*exp(x) )
    elseif 0.1 <= x < 0.5
        exp(w0) / ( exp(w0) + multi*0.5*exp(x) )        
    else
        exp(w0) / ( exp(w0) + multi*exp(x) )       
    end        
end

num_test_statistics = length(test_statistics)
posterior = Array{Float64}(undef, num_test_statistics)

for i in 1:num_test_statistics
    ts = test_statistics[i]
    prior_val = prior_fn(prior_list[i])
    null_val = null_distr
    mix_val = mixture_pdf(ts)

    # if mixture distr equals 0, then just return a 0 posterior
    if mix_val == zero(mix_val)
        posterior[i] = 0.0
        continue
    end

    fdr = null_val / mix_val
    p1 = 1 - prior_val * fdr
    posterior[i] = p1
end

#----------------------------------Formatting time-------------------------------------------
#List where the edges will be stored
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
    print(weight)
    print('\n')
    if weight < p_val
        break
    end
    edge1 = parse(Int,edge1[2:end]) + 1
    edge2 = parse(Int,edge2[2:end]) + 1
    matrix[edge1,edge2] = 1
    matrix[edge2,edge1] = 1 
end

#writes the file to a matrix
num_cells = Values[4]
# type = Values[5]
name = join(["Eb_matrix_", num_cells, ".csv"])
writedlm( name,  matrix, ',')
