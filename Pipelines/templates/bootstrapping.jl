#Here we will use bootstrapping for the Empirical Bayes method

# include("../../julia_environment/src/julia_environment.jl")

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

#Data_matrix = The input matrix
#sample_num = Number of the interation
#prior_matrix = Prior matrix generated
#p_val = p_value used
function EB_running(data_matrix,sample_num,prior_matrix,p_val)
    num_col = size(data_matrix)[2] #Get number of rows
    samples = rand(1:num_col,2000) #Get a sample of the columns
    bs_mat = data_matrix[:,samples] #Bootstrapped matrix
    #Now the tricky bit, making a matrix it will be happy with
    i = sample_num
    bs_name = join(["sample_", i, ".txt"])
    open(joinpath([path,bs_name]), "w") do file
        println(file,"0\t0\t0\t50\t50\t50\t100\t100\t100\t500\t600\t500\t\n")
        writedlm(file, [gene_names bs_mat])
    end

    #Network inference time
    genes = get_nodes(joinpath([path,bs_name]));
    inference = PUCNetworkInference()
    network_EB = InferredNetwork(inference, genes);

    #Edge list
    edge_list = network_EB.edges
    test_statistics = [e.weight for e in edge_list]

    #Get the prior lists
    prior_list = [ get(prior_matrix, to_index(e.nodes), 0) for e in edge_list ]

    #Running EB
    num_bins = 5
    distr = :Normal
    proportion_to_keep = 0.9
    tail = :two
    w0 = 2.2
    posteriors = empirical_bayes(test_statistics, prior_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)

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

    p_val = p_val

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
    return matrix
end

Values = ARGS

data = Values[1] #Read in data
data_matrix = readdlm(data, '\t', Float64) #Put into matrix
#Gets the gene names
gene_names = []
num_genes = size(data_matrix)[1]
for i in 0:(num_genes-1)
    cur_gene = string("T", i)
    push!(gene_names,cur_gene)
end

#Prior information
priors = Values[2]
priors = CSV.read(priors, DataFrame, header=false)
#We convert this into a matrix
prior_matrix = Dict()
for i in 1:nrow(priors)
    edge1 = parse(Int,priors[i,1][2:end]) + 1
    edge1_name = join(["T", string(edge1)])
    edge2 = parse(Int,priors[i,2][2:end]) + 1
    edge2_name = join(["T", string(edge2)])
    prior_matrix[(edge1_name,edge2_name)] = priors[i,3]
end

#Creating a temporary folder
path = "temp_folder"
mkdir(path)

all_values = Matrix{Int64}[]
N = 4
p_val = parse(Float64,Values[3])
Threads.@threads for i = 1:N   
    push!(all_values,EB_running(data_matrix,i,prior_matrix,0.9))
end

id = Values[4]
matrix = sum(all_values)
name = join(["EB_", id, "_bootstrapping.csv"])
writedlm(name .> 2/N,  matrix, ',')

# rm(path,recursive=true)

