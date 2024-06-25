#Here is where I will write the functions which will be tested with the Empirical Bayes script

#---MODULES USED---
using DataFrames, CSV
using NetworkInference
using EmpiricalBayes
using EvalMetrics
using Distributions

#---END OF MODULES---

#A very simple test to determine if I can run what I want to run
function basic_test(x)
    return (x+1)
end

#Some functions needed for Empirical Bayes to run
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

#Now we can test the prior function
function get_priors(data, prior_data)
   #Gets the network
   genes = get_nodes(data);

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

   return priors
end

#Defines the distribution
function get_distr(d::Symbol)
    distrs = Dict(
        :Gamma => Gamma,
        :Normal => Normal
    )
    return distrs[d]
end

#This function runs Empirical Bayes
function run_EB(data,prior_data,w0,proportion_to_keep,multi)
    genes = get_nodes(data);
    num_cells = ncol(DataFrame(CSV.File(data, header= true))) - 1
    if num_cells < 2000
        inference = MINetworkInference()
        distr = :Gamma
    else
        inference = PUCNetworkInference()
        distr = :Normal
    end
    network_EB = InferredNetwork(inference, genes);
    edge_list = network_EB.edges
    test_statistics = [e.weight for e in edge_list]

    priors = get_priors(data,prior_data)
    priors = [ get(priors, to_index(e.nodes), 0) for e in edge_list ]
    num_bins = 5

    #Runs the EB stuff
    midpoints, counts, bin_width = discretize_test_statistics(test_statistics, num_bins)
    null_distr = fit_null_distribution(midpoints, counts, num_bins, bin_width, proportion_to_keep, get_distr(distr))
    mixture_pdf = fit_mixture_distribution(midpoints, counts, bin_width)

    #From here we will be analysing the Posterior calculation
    null_pdf(x) = pdf(null_distr, x)

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
    
    #Here is where the magic happens
    for i in 1:num_test_statistics
        ts = test_statistics[i]
        prior_val = prior_fn(priors[i])
        null_val = null_pdf(ts)
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
    #List where the edges will be stored
    eb_edges = Array{Edge}(undef, length(edge_list))
    for i in 1:length(edge_list)
        nodes = edge_list[i].nodes
        eb_edges[i] = Edge(nodes, posterior[i])
    end

    p_val = 0.95

    #Here is the final touch
    num_genes = length(genes)
    matrix = zeros(Int,num_genes,num_genes)
    for i = 1:length(eb_edges)
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
    return matrix
end

#Let's now calculate the AUPR based on the returned function:
function scores(truth,matrix)
    guess = convert(Matrix{Float64},matrix)
    actual = Matrix(truth)
    return au_prcurve(vec(actual),vec(guess))
end

#We now test out the function
data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/Michael/EB_formatted/formatted_data.txt"
prior = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/Michael/EB_formatted/ground_truth1-50_2.2.csv"
w0 = 2.2
to_keep = 0.8
multi = 10
# matrix = run_EB(data,prior,w0,to_keep,multi)

#Function to assess EvalMetrics
function metrics(data,prior,w0,to_keep,multi,truth)
    matrix = run_EB(data,prior,w0,to_keep,multi)
    guess  = convert(Matrix{Float64},matrix)
    actual = DataFrame(CSV.File(truth, header= false))
    actual = Matrix(actual)
    return au_prcurve(vec(actual),vec(guess))
end

metrics(data,prior,w0,to_keep,multi,prior)