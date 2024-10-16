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
full_data = Values[2]
id = Values[5]
zero_data = Values[3]
genie_data = Values[4]
nlnet_data = Values[5]
w0 = 2.2
multi = 10
p_val = parse(Float64,Values[6])

genes = get_nodes(data);
network_EB = InferredNetwork(inference, genes);
#-------Here I am just getting the PIDC network-----------------------
adjacency_matrix, labels_to_ids, ids_to_labels = get_adjacency_matrix(network_EB, 0.1)
name = join(["PIDC_", id,".csv"])
writedlm( name,  Int.(adjacency_matrix), ',')

edge_list = network_EB.edges;
# #This here will spit out the test statistics needed
test_statistics = [e.weight for e in edge_list]

#------------------------------------------Priors----------------------------------------------

full_file = DataFrame(CSV.File(full_data, header= false))
full_file = Matrix(full_file)
zero_file = DataFrame(CSV.File(zero_data, header= false))
zero_file = Matrix(zero_file)
genie_file = DataFrame(CSV.File(genie_data, header= false))
genie_file = Matrix(genie_file)
# nlnet_file = DataFrame(CSV.File(nlnet_data, header= false))
# nlnet_file = Matrix(nlnet_file)


#Now generates a dictionary of the prior information
num_genes = length(genes)
# Generate the names of the gene
gene_names = []
for i in 0:(num_genes-1)
    cur_gene = string("T", i)
    push!(gene_names,cur_gene)
end
#Builds the dictionary
full_priors = Dict()
zero_priors = Dict()
genie_priors = Dict()
nlnet_priors = Dict()
for i in gene_names
    for j in gene_names
        edge1 = parse(Int,i[2:end]) + 1
        edge2 = parse(Int,j[2:end]) + 1
        full_priors[(i,j)] = full_file[edge1,edge2]*2.2
        zero_priors[(i,j)] = zero_file[edge1,edge2]*2.2
        genie_priors[(i,j)] = genie_file[edge1,edge2]*2.2
        # nlnet_priors[(i,j)] = nlnet_file[edge1,edge2]*2.2
    end
end

#We convert this into a matrix
# full_matrix = Dict()
# for i in 1:nrow(prior_file)
#     edge1 = parse(Int,prior_file[i,1][2:end]) + 1
#     edge1_name = join(["T", string(edge1)])
#     edge2 = parse(Int,prior_file[i,2][2:end]) + 1
#     edge2_name = join(["T", string(edge2)])
#     full_matrix[(edge1_name,edge2_name)] = prior_file[i,3]
# end

full_list  = [ get(full_priors, to_index(e.nodes), 0) for e in edge_list ]
genie_list = [ get(genie_priors, to_index(e.nodes), 0) for e in edge_list ]
zero_list  = [ get(zero_priors, to_index(e.nodes), 0) for e in edge_list ]
# nlnet_list = [ get(nlnet_priors, to_index(e.nodes), 0) for e in edge_list ]


#----------------------------------Where the magic happens--------------------------------------
midpoints, counts, bin_width = discretize_test_statistics(test_statistics, num_bins)
null_distr = 0.5 #Uniform distribution
mixture_pdf = fit_mixture_distribution(midpoints, counts, bin_width) #Mixture is important


#--------------------------------Method to test now-------------------------------------------

# Prior function
multi = 10
function prior_fn(x) 
    if x == 0
        exp(w0) / ( exp(w0) + exp(x) )
    elseif 0 < x < 0.22
        exp(w0) / ( exp(w0) + 0*exp(x) )
    elseif 0.22 <= x < 1.1
        exp(w0) / ( exp(w0) + multi*0.5*exp(x) )        
    else
        exp(w0) / ( exp(w0) + multi*exp(x) )       
    end        
end

# prior_fn(x) = exp(w0) / ( exp(w0) + exp(x) )

num_test_statistics = length(test_statistics)
full_posteriors = Array{Float64}(undef, num_test_statistics)
zero_posteriors = Array{Float64}(undef, num_test_statistics)
genie_posteriors = Array{Float64}(undef, num_test_statistics)
nlnet_posteriors = Array{Float64}(undef, num_test_statistics)

for i in 1:num_test_statistics
    ts = test_statistics[i]
    full_val = prior_fn(full_list[i])
    genie_val = prior_fn(genie_list[i])
    zero_val = prior_fn(zero_list[i])
    # nlnet_val = prior_fn(nlnet_list[i])
    null_val = null_distr
    mix_val = mixture_pdf(ts)

    # if mixture distr equals 0, then just return a 0 posterior
    if mix_val == zero(mix_val)
        full_posteriors[i] = 0.0
        zero_posteriors[i] = 0.0
        genie_posteriors[i] = 0.0
        nlnet_posteriors[i] = 0.0
        continue
    end

    fdr = null_val / mix_val
    full_posteriors[i] = 1 - full_val * fdr
    genie_posteriors[i] = 1 - genie_val * fdr
    zero_posteriors[i] = 1 - zero_val * fdr
    # nlnet_posteriors[i] = 1 - nlnet_val*fdr
end

#----------------------------------Formatting time-------------------------------------------
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
name = join(["Eb_fullmatrix_", id, ".csv"])
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
name = join(["Eb_geniematrix_", id, ".csv"])
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
name = join(["Eb_zeromatrix_", id,".csv"])
writedlm( name,  matrix, ',')

#------NLNET--------
# #List to hold edges
# edges = Array{Edge}(undef, length(edge_list))
# #We also want to get weighted values
# weighted = zeros(length(edge_list))
# #Gets the edges
# for i in 1:length(edge_list)
#     nodes = edge_list[i].nodes

#     edges[i] = Edge(nodes, nlnet_posteriors[i])
#     weighted[i] = edges[i].weight
# end

# permvec = sortperm(weighted)
# #We now want to create a matrix which holds all the values
# matrix = zeros(Int,num_genes,num_genes)
# for i in reverse(permvec)
#     edge1 = edges[i].nodes[1].label
#     edge2 = edges[i].nodes[2].label
#     weight = edges[i].weight
#     if weight < p_val
#         break
#     end
#     edge1 = parse(Int,edge1[2:end]) + 1
#     edge2 = parse(Int,edge2[2:end]) + 1
#     matrix[edge1,edge2] = 1
#     matrix[edge2,edge1] = 1 
# end
# #writes the file to a matrix
# name = join(["Eb_nlnetmatrix_", id,".csv"])
# writedlm( name,  matrix, ',')
