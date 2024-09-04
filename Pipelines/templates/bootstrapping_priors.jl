#Here we will use bootstrapping for the Empirical Bayes method

# include("../../julia_environment/src/julia_environment.jl")

using NetworkInference
using EmpiricalBayes
using DelimitedFiles
using DataFrames, CSV
using Distributions

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

function get_distr(d::Symbol)
    distrs = Dict(
        :Gamma => Gamma,
        :Normal => Normal
    )
    return distrs[d]
end

#Data_matrix = The input matrix
#sample_num = Number of the interation
#prior_matrix = Prior matrix generated
#p_val = p_value used
function EB_running(data_matrix,sample_num,prior_matrix,zero_matrix,p_val)
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
    zero_list = [ get(zero_matrix, to_index(e.nodes), 0) for e in edge_list ]

    #Running EB
    num_bins = 5
    distr = :Normal
    w0 = 2.2


    #-------------------------------Running EB the way that is needed-----------------------------------------------
    # posteriors = empirical_bayes(test_statistics, prior_list, num_bins, distr, proportion_to_keep = proportion_to_keep, tail = tail, w0 = w0)

    midpoints, counts, bin_width = discretize_test_statistics(test_statistics, num_bins)
    #MLE estimates
    norm_mean = (1/length(test_statistics))*sum(test_statistics)
    norm_sigma = sqrt( (1/length(test_statistics))* sum((test_statistics .- norm_mean).^2) )
    distr = get_distr(distr)
    null_distr = distr(norm_mean,norm_sigma)
    mixture_pdf = fit_mixture_distribution(midpoints, counts, bin_width)


    #From here we will be analysing the Posterior calculation
    null_pdf(x) = pdf(null_distr, x)
    multi = 10

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
    posteriors_prior = Array{Float64}(undef, num_test_statistics)
    posteriors_zero = Array{Float64}(undef, num_test_statistics)

    #Here is where the magic happens
    for i in 1:num_test_statistics
        ts = test_statistics[i]
        prior_val = prior_fn(prior_list[i])
        zero_val = prior_fn(zero_list[i])
        null_val = null_pdf(ts)
        mix_val = mixture_pdf(ts)

        # if mixture distr equals 0, then just return a 0 posterior
        if mix_val == zero(mix_val)
            posteriors_prior[i] = 0.0
            posteriors_zero[i] = 0.0
            continue
        end

        fdr = null_val / mix_val
        p1 = 1 - prior_val * fdr
        posteriors_prior[i] = p1
        posteriors_zero[i] = 1 - zero_val*fdr
    end
    #-------------------------------Finished running EB-------------------------------------------------------------


    #List to hold edges
    edges = Array{Edge}(undef, length(edge_list))
    #We also want to get weighted values
    weighted = zeros(length(edge_list))
    #Gets the edges
    for i in 1:length(edge_list)
        nodes = edge_list[i].nodes

        edges[i] = Edge(nodes, posteriors_prior[i])
        weighted[i] = edges[i].weight
    end

    p_val = p_val

    permvec = sortperm(weighted)
    matrix_priors = zeros(Int,num_genes,num_genes)
    for i in reverse(permvec)
        edge1 = edges[i].nodes[1].label
        edge2 = edges[i].nodes[2].label
        weight = edges[i].weight
        if weight < p_val
            break
        end
        edge1 = parse(Int,edge1[2:end]) + 1
        edge2 = parse(Int,edge2[2:end]) + 1
        matrix_priors[edge1,edge2] = 1
        matrix_priors[edge2,edge1] = 1 
    end

    #--------------Now we get the zero values---------------
    #List to hold edges
    edges = Array{Edge}(undef, length(edge_list))
    #We also want to get weighted values
    weighted = zeros(length(edge_list))
    #Gets the edges
    for i in 1:length(edge_list)
        nodes = edge_list[i].nodes

        edges[i] = Edge(nodes, posteriors_zero[i])
        weighted[i] = edges[i].weight
    end

    p_val = p_val

    permvec = sortperm(weighted)
    matrix_zero = zeros(Int,num_genes,num_genes)
    for i in reverse(permvec)
        edge1 = edges[i].nodes[1].label
        edge2 = edges[i].nodes[2].label
        weight = edges[i].weight
        if weight < p_val
            break
        end
        edge1 = parse(Int,edge1[2:end]) + 1
        edge2 = parse(Int,edge2[2:end]) + 1
        matrix_zero[edge1,edge2] = 1
        matrix_zero[edge2,edge1] = 1 
    end

    return matrix_priors, matrix_zero
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
zero_matrix = Dict()
for i in 1:nrow(priors)
    edge1 = parse(Int,priors[i,1][2:end]) + 1
    edge1_name = join(["T", string(edge1)])
    edge2 = parse(Int,priors[i,2][2:end]) + 1
    edge2_name = join(["T", string(edge2)])
    prior_matrix[(edge1_name,edge2_name)] = priors[i,3]*2.2
    zero_matrix[(edge1_name,edge2_name)] = 0
end

#Creating a temporary folder
path = "temp_folder"
mkdir(path)

prior_values = Matrix{Int64}[]
zero_values = Matrix{Int64}[]
N = 50
p_val = parse(Float64,Values[3])
Threads.@threads for i = 1:N  
    matrix_priors, matrix_zeros = EB_running(data_matrix,i,prior_matrix,zero_matrix,0.9)
    push!(prior_values,matrix_priors)
    push!(zero_values,matrix_zeros)
end

id = Values[4]
matrix_prior = sum(prior_values)
name = join(["EB_", id, "_bootstrapping.csv"])
writedlm( name, Int64.(matrix_prior .> 25/N), ',')

matrix_zero = sum(zero_values)
name = join(["EB_", id, "_bootstrapping_zero.csv"])
writedlm( name, Int64.(matrix_zero .> 25/N), ',')

# rm(path,recursive=true)

