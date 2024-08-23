// This file will be running GENIE3 on all the datasets that are needed
//We aim to perform this in a six step process

//The scripts we will use
genie_research = file( params.genie_research )
genie_script   = file( params.genie_script )
genie_con_gene = file( params.genie_con_gene )
eb_run         = file( params.Emp_bay )
input_fix      = file( params.input_fix )
metrics_mich   = file( params.metrics_mich )

//A tuple of values which will be used for the Empirical Bayes run
eb_tuple = tuple(0.75, 0.80)

workflow{

    //STEP 1: Get all data needed into the appropriate tuples

    //We will get a tuple with the following structure:
    //(ID, training dataset, testing dataset, ground truth network)

    //Ground truth networks
    def groundtruth = Channel.fromPath( "../Actual_data_tests/Research_data/inputs/Curated/HSC/HSC-2000-*/refNetwork.csv" )
    //The training data
    def training_expressions = Channel.fromPath( "../Actual_data_tests/Research_data/inputs/Curated/HSC/HSC-2000-*/ExpressionData.csv" )
    //The testing data
    def testing_expressions = Channel.fromPath( "../Actual_data_tests/Research_data/inputs/Curated/HSC/HSC-2000-*-50/ExpressionData.csv" )
    //The values of pi0 we will be testing for 2.2 = 0.9, 2.944 = 0.95
    def pi0 = Channel.of(2.2,2.944)    
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    training_expressions
        .map { [it.toString().split('2000-')[1].split('/ExpressionData.csv')[0], it] }
        .set { expressiontrain_key }
    testing_expressions
        .map { [it.toString().split('2000-')[1].split('-50/ExpressionData.csv')[0], it] }
        .set { expressiontest_key }
    groundtruth
        .map { [it.toString().split('2000-')[1].split('/refNetwork.csv')[0], it] }
        .set { groundtruth_key }
    ch_inputs = expressiontrain_key
        .join(expressiontest_key, by:0)
        .join(groundtruth_key, by:0)


    //STEP 2: Run all clean ups on it (including getting reference networks)

    //This creates a tuple with the following values:
    //(ID, converted training data, converted testing data, truth matrix, gene names)
    convert_ch = CONVERT(
        genie_research,
        ch_inputs
    )
    ch_inputs.view()

    //STEP 3: Run EB on the 50 file (using no priors)
    
    //This will create a tuple with the following values
    // (ID, converted testing data, truth matrix, gene names, genie3 outputs)
    // genie3_ch = GENIE(
    //     genie_script,
    //     convert_ch
    // )

    // //STEP 4: Normalise the output and put reference network in a matrix

    // //This will create a tuple with the following values
    // // (ID, converted testing data, truth matrix, gene names, prior matrix)

    // prior_ch = GENIE_CONVERSION(
    //     genie_con_gene,
    //     genie3_ch
    // )


    // //STEP 5: With the reference network, run EB again, this time on the 70 data using 50 prior

    // //Firstly convert the training data into something worthwhile
    // prior_ch = EB_CONVERT(
    //     input_fix,
    //     prior_ch
    // )

    // //Secondly combine with the pi0 values
    // prior_ch = prior_ch.combine(pi0)

    // //Run Empirical Bayes
    // //This will return a tuple of the following nature
    // // (Id, true network, gene list, pi0 value used, the EB network)
    // eb_ch = EB_RUN(
    //     eb_run,
    //     eb_tuple,
    //     prior_ch
    // )

    // //STEP 6: Return all of the metrics
    // //Returns 18 metrics
    // METRICS(
    //     metrics_mich,
    //     eb_ch
    // )

}

//Here we convert the given reference network into something useful
//It will also help convert the data files into something workable
//It also gets a list of genes which will be used in later calculations
process CONVERT{

    publishDir "${params.outdir}/prior_genie/Converted"

    input:
    path script 
    tuple val(id), path(training, stageAs: 'training.csv'), path(testing, stageAs: 'testing.csv'), path(truth)

    output:
    tuple val(id), path("training${id}.txt"), path("testing${id}.txt"),
        path("truth${id}.csv"), path("genes${id}.txt"), emit: tuple_converter

    """
    python3 ${script} ${training} ${testing} ${truth} ${id}
    """
}

//Here we run GENIE3 on the testing data
//We do no normalising here
process GENIE{

    cpus 4
    memory '4 GB'
    container 'genie3:latest'
    publishDir "${params.outdir}/prior_genie/genie_outputs"

    input:
    path genie_script
    tuple val(id), path(training), path(testing), path(truth), path(genes)

    output:
    tuple val(id), path(testing), path(truth), path(genes), path("genie${id}.txt")

    script:

    """
    Rscript ${genie_script} ${training} ${id} > genie${id}.txt
    """

}

//Here we will deal with formatting of the output
//This will normalise the data into something which can be used for the output
process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/prior_genie/priors", mode: 'copy'

    input:
    path script
    tuple val(id), path(testing), path(truth), path(genes), path(genie_output)

    output:
    tuple val(id), path(testing), path(truth), path(genes), path("matrix_genie_normalised${id}.csv")

    script:
    """
    python3 ${script} ${genie_output} ${genes} ${id}
    """
}

//Here we convert the file into something useable
process EB_CONVERT{

    publishDir "${params.outdir}/prior_genie/EB_formatted"

    input:
    path script
    tuple val(id), path(testing), path(truth), path(genes), path(prior)

    output:
    tuple val(id), path('formatted_data.txt'), path(truth), 
        path(genes), path(prior), emit: tuple_converter

    script:

    """
    python3 ${script} ${testing} > formatted_data.txt
    """    

}

//Here we run Empirical Bayes
process EB_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/prior_genie/EB_output"

    input:
    path script
    tuple val(p_val), val(keep)
    tuple val(id), path(testing), path(truth), path(genes), path(prior), val(pi0)

    output:
    tuple val(id), path(truth), path(genes), val(pi0), path("Eb_matrix_${id}_${pi0}.csv")

    script:

    """
    julia ${script} ${testing} ${p_val} ${prior} ${id} ${keep} ${pi0}
    """
}

//Runs the metrics
process METRICS{

    container 'metrics:latest'
    publishDir "${params.outdir}/prior_genie/metrics"

    input:
    path script 
    tuple val(id), path(truth), path(genes), val(pi0), path(guess)

    output:
    path "metrics_${id}_${pi0}.txt"

    """
    python3 ${script} ${guess} ${truth} ${genes} ${pi0} ${id}
    """
}