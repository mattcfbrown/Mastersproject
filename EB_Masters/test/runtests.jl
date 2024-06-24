#Here I will write some Julia testing code for Empirical Bayes
#Step 1 is working out how to get testing to work

using Test
using DelimitedFiles

#We are firstly aiming to see if we can run tests
@testitem "Can I write tests?" begin
    x = "Hello world"
    @test length(x) == 11
end

#Now we can test out some of the EB packages we need to run
@testitem "Can I write tests in another function" begin
    x = 1
    @test basic_test(x) == 2
end

#Let's firstly test to make sure we are inputting prior information correctly
@testitem "Prior tests" begin
    data = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/Michael/EB_formatted/formatted_data.txt"
    prior = "/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Prior_testing/results/Michael/EB_formatted/ground_truth1-50_2.2.csv"
    @test get_priors(data,prior)[("T0", "T2")] == 2.2 #Goes forwards
    @test get_priors(data,prior)[("T2", "T0")] == 2.2 #Goes backwards
    @test get_priors(data,prior)[("T8", "T2")] == 0.0 #zero
    @test get_priors(data,prior)[("T11", "T4")] == 0.0 #Another zero
    @test get_priors(data,prior)[("T15", "T15")] == 2.2 #Diagonal
end

#I want to run the following tests:

#Good prior
#Bad prior

#For each have three thresholds
#It is better than garbage
#It is better than average
#It is better than Good

#This means I will need to write out six new test cases
