#Packages used
using NetworkInference
using Plots
using EmpiricalBayes
using Distributions

#Data used
# data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/PIDC_testing/Data/shuffled.txt"
# data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/PIDC_testing/Data/indepenent.txt"
data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/thesis/EB_formatted/formatted_data_250cells.txt"
inference = PIDCNetworkInference()

# #This file will aim to do the following:

# #Step 1: Run Network Inference using PIDC method
genes = get_nodes(data);
network_EB = InferredNetwork(inference, genes);
edge_list = network_EB.edges;
# #This here will spit out the test statistics needed
test_statistics = [e.weight for e in edge_list]

# #We remove the top 10% of values
threshold_index = convert(Int64,round(length(test_statistics)*0.1))
threshold = test_statistics[threshold_index]
ten_removed = test_statistics[test_statistics .< threshold]

#Step 2: Get a normal null dsitribution estimate based off the values
#This method comes from EmpiricalBayes Package
#I have already removed the top 10% of statistics, and as such, will not need to further improve on it 
#Firstly some variables needed for below
num_bins = 10 #Number of bins, reflects the histogram above
# proportion_to_keep = 1.0 #We keep it all as we did it before

#For the distributions
function get_distr(d::Symbol)
    distrs = Dict(
        :Gamma => Gamma,
        :Normal => Normal
    )
    return distrs[d]
end

#Do part 1
midpoints, counts, bin_width = discretize_test_statistics(ten_removed, num_bins)
#Now part 2
# null_distr_norm = fit_null_distribution(midpoints, counts, num_bins, bin_width, proportion_to_keep, Normal)
# null_distr_gamm = fit_null_distribution(midpoints, counts, num_bins, bin_width, proportion_to_keep, Gamma)

#This ius a function we can use for plotting
# null_norm(x) = pdf(null_distr_norm, x)
# null_gamm(x) = pdf(null_distr_gamm, x)

#Step 3: Get histogram of the values used, plot the normal over the top
num_edges = length(test_statistics)
b_range = range(0,2, length = 30)

# plot(sort(test_statistics), (1:num_edges)./num_edges,
# title = "PIDC", label = "PIDC CDF")
# plot!(range(0,2, length = 50),
# range(0,1, length = 50),
# label = "Uniform CDF")
# xlabel!("c")
# ylabel!("F(c)")

histogram(test_statistics, normalize=:pdf, label = "PIDC", color = :blue, bins=b_range, 
title = "PIDC pdf")
plot!([0,2],[0.5,0.5],
label = "Uniform")
xlabel!("c")
ylabel!("f(c)")
# plot!(null_norm, label = "Normal", color = :red)
# plot!(null_gamm, label = "Gamma", color = :blue)

#-----------------------------------------------------------------
#A new method to be developed here

#Step 1: Get the PUC scores:
# inference = PUCNetworkInference()
# genes = get_nodes(data);
# network_EB = InferredNetwork(inference, genes);
# #This here gives us two labels and it's corresponding weight
# # print(network_EB.edges[1].nodes[1].label)
# # print(network_EB.edges[1].nodes[2].label)
# # print(network_EB.edges[1].weight)
# Weights = Dict()
# num_edges = length([e.weight for e in network_EB.edges])
# for i in 1:num_edges
#     Node1 = network_EB.edges[i].nodes[1].label
#     Node2 = network_EB.edges[i].nodes[2].label
#     value = network_EB.edges[i].weight

#     if haskey(Weights,Node1)
#         push!(Weights[Node1],value)
#     else
#         merge!(Weights, Dict(Node1 => [value]))
#     end

#     if haskey(Weights,Node2)
#         push!(Weights[Node2],value)
#     else
#         merge!(Weights, Dict(Node2 => [value]))
#     end
# end

# #Step 2: Get the pdf score for each situation
# #Holds the eventual test statistics
# test_stats = zeros(num_edges)
# for i in 1:2
#     Node1 = network_EB.edges[i].nodes[1].label
#     Node2 = network_EB.edges[i].nodes[2].label
#     value = network_EB.edges[i].weight
#     test = pdf( fit(Normal,Weights[Node1]), value)
#     print(test)
# end
