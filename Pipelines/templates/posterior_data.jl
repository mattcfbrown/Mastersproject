#This file aims to output two realisations
#A density graph for the data + Posterior
#A density graph for their division, with pi0 values where values would accepted

#Firstly let's run Empirical Bayes
using NetworkInference
using Distributions
using EmpiricalBayes
using DelimitedFiles
using DataFrames, CSV
using Plots

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

# Values = ARGS

#All the values needed to perform the calculations
#This here gets the test statistics
data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Empirical_Bayes_testing/results/Information_Measures/formatted_data.txt"
genes = get_nodes(data);
num_cells = ncol(DataFrame(CSV.File(data, header= true))) - 1
if num_cells < 500
    inference = MINetworkInference()
    distr = :Gamma
else
    inference = PUCNetworkInference()
    distr = :Normal
end
network_EB = InferredNetwork(inference, genes);
edge_list = network_EB.edges
test_statistics = [e.weight for e in edge_list]


num_bins = 5
proportion_to_keep = 0.9                                   #Change this
p_val = 0.9                                                #Change this
tail = :two
w0=2.2
null_value = -Inf

function get_distr(d::Symbol)
    distrs = Dict(
        :Gamma => Gamma,
        :Normal => Normal
    )
    return distrs[d]
end

#This here is what is going to help us
midpoints, counts, bin_width = discretize_test_statistics(test_statistics, num_bins)
null_distr = fit_null_distribution(midpoints, counts, num_bins, bin_width, proportion_to_keep, get_distr(distr))
mixture_pdf = fit_mixture_distribution(midpoints, counts, bin_width)

#We will now work using this as a basis
null_pdf(x) = pdf(null_distr, x)


num_test_statistics = length(test_statistics)
#We hold the values here
null_val = Array{Float64}(undef, num_test_statistics)
mix_val = Array{Float64}(undef, num_test_statistics)
total = Array{Float64}(undef, num_test_statistics)

#This is an edited for loop to hold our data
for i in 1:num_test_statistics
    ts = test_statistics[i]
    null_val[i] = null_pdf(ts)
    mix_val[i] = mixture_pdf(ts)
    total[i] = 1 - null_val[i]/mix_val[i]
end

#Here we get the posterior values which are the thresholds
alpha = [0.5,0.8,0.9,0.95]
values = 1.0.-(1-p_val)./alpha

#Plot the graphs
b_range = range(-0.05, 0.1, length=50)
histogram([null_val mix_val],                                               #The data, null and mixture distributions
            label = ["Null dsitribution" "Mixture Distribution"],           #The labels used
            opacity=[0.6 0.6],                                              #Makes the data see thrrough
            normalize=:pdf,                                                 #Normalises it
            bins = b_range)                                                 #Create bins
xlims!(-0.1, 0.1)
title!("Distribution of the likelihood (null) and data (denominator)",
        titlefont = font(12,"Computer Modern"))
xlabel!("x")
ylabel!("density")

savefig("plot1.pdf")


#Now do it for the division
histogram(total,                                                            #The data
            label = "Null/mixture",                                         #The labels used
            normalize=:pdf)                                                 #Normalises it    
title!("Distribution of the likelihood/mixture",
    titlefont = font(12,"Computer Modern"))
xlabel!("x")
ylabel!("density")
vline!([values[1]], color = [:red], label =string(alpha[1]))
vline!([values[2]], color = [:green], label =string(alpha[2]))
vline!([values[3]], color = [:magenta], label =string(alpha[3]))
vline!([values[4]], color = [:tan], label =string(alpha[4]))

savefig("plot2.pdf")