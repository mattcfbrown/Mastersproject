using NetworkInference
using EmpiricalBayes

Values = ARGS


#Gets the network
genes = get_nodes(Values[1]);
network = InferredNetwork(PUCNetworkInference(), genes);

#Perform EmpiricalBayes on it 
eb_network = empirical_bayes(network, 5, :Gamma, tail = :two, w0 = 2.2)
println("Here are the weightings\n")
print([e.weight for e in eb_network.edges])

