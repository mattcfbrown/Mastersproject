// This program here will run four different tests:
//It will run no priors and full priors

//Scripts used
research_con  = file( params.research_con )
input_fix     = file( params.input_fix )
eb_run        = file( params.Emp_bay )
metrics_mich  = file( params.metrics_mich )

//Some of the files I will use (as so I can have a toy example):
toy_truth     = file( params.toy_truth )
toy_data      = file( params.toy_data )
zero          = file( params.zero )


//A tuple of values which will be used for the Empirical Bayes run
eb_tuple = tuple(0.75, 0.80)


//Main program here:
workflow{

    //STEP 1: Get all data needed into the appropriate tuples

    //We will get a tuple with the following structure:
    //(ID, dataset, ground truth network, pi0)

    //This is what I will use for all of the data if needed
    //This holds all the information:
    def groundtruth = Channel.fromPath( "../Actual_data_tests/Research_data/inputs/Curated/HSC/HSC-2000-?-50/refNetwork.csv" )
    def expressions = Channel.fromPath( "../Actual_data_tests/Research_data/inputs/Curated/HSC/HSC-2000-?-50/ExpressionData.csv" )
    def pi0 = Channel.of(2.2,2.944)
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    expressions
        .map { [it.toString().split('2000-')[1].split('/ExpressionData.csv')[0], it] }
        .set { expression_key }
    groundtruth
        .map { [it.toString().split('2000-')[1].split('/refNetwork.csv')[0], it] }
        .set { groundtruth_key }
    ch_inputs = expression_key
        .join(groundtruth_key, by:0)
        .combine(pi0)

    //STEP 2: Run all clean ups on it (including getting reference networks)

    //This creates a tuple with the following values:
    //(ID, converted data, truth matrix, gene names, pi0)

    output_ch = CONVERT(
        research_con,
        ch_inputs
    )

    //STEP 3: With the reference network, run EB, Choosing the type of prior information you wish to use

    //Firstly convert the training data into something worthwhile
    formatted_data = EB_CONVERT(
        input_fix,
        output_ch
    )

    // We can now add the prior information
    full_ch  = formatted_data.combine(Channel.of('full'))
    empty_ch = formatted_data.combine(Channel.fromPath(zero)).combine(Channel.of('zero'))

    //Run Empirical_Bayes for full priors:
    full_ch = EB_RUN_FULL(
        eb_run,
        eb_tuple,
        full_ch     
    )

    //Run EB on no priors
    empty_ch = EB_RUN_EMPTY(
        eb_run,
        eb_tuple,
        empty_ch   
    )

    //We get two channels of the following types:
    // (ID, ground truth, gene list, pi0 value, matrix, type of prior)

    //STEP 4: Return all of the metrics
    //Returns 36 metrics

    //Now run the metrics (Full)
    METRICS_FULL(
        metrics_mich,
        full_ch
    )

    // //Empty
    METRICS_NULL(
        metrics_mich,
        empty_ch        
    )
}


//Converts the file into something workable
process CONVERT{

    publishDir "${params.outdir}/Michael/Converted"

    input:
    path script 
    tuple val(id), path(data), path(truth), val(pi0)

    output:
    tuple val(id), path("expression${id}_${pi0}.txt"), path("ground_truth${id}_${pi0}.csv"), 
        path("genes${id}_${pi0}.txt"), val(pi0), emit: tuple_converter

    """
    python3 ${script} ${data} ${truth} ${id} ${pi0}
    """
}

//Converts the Empirical Bayes into something workable
process EB_CONVERT{

    publishDir "${params.outdir}/Michael/EB_formatted"

    input:
    path script
    tuple val(id), path(data), path(truth), path(genes), val(pi0)

    output:
    tuple val(id), path('formatted_data.txt'), path(truth), 
        path(genes), val(pi0), emit: tuple_converter

    script:

    """
    python3 ${script} ${data} > formatted_data.txt
    """    

}

//Runs Empirical Bayes for full priors
process EB_RUN_FULL {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/Michael/EB_output"

    input:
    path script
    tuple val(p_val), val(keep)
    tuple val(id), path(data), path(truth), path(genes), val(pi0), val(type)

    output:
    tuple val(id), path(truth), path(genes), val(pi0), path("Eb_matrix_${id}_${pi0}_${type}.csv"), val(type)

    script:

    """
    julia ${script} ${data} ${p_val} ${truth} ${id} ${keep} ${pi0} ${type}
    """
}

//Runs EB for no priors
process EB_RUN_EMPTY {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/Michael/EB_output"

    input:
    path script
    tuple val(p_val), val(keep)
    tuple val(id), path(data), path(truth), path(genes), val(pi0), path(prior), val(type)

    output:
    tuple val(id), path(truth), path(genes), val(pi0), path("Eb_matrix_${id}_${pi0}_${type}.csv"), val(type)

    script:

    """
    julia ${script} ${data} ${p_val} ${prior} ${id} ${keep} ${pi0} ${type}
    """
}

//Metrics for the full version
process METRICS_FULL{

    container 'metrics:latest'
    publishDir "${params.outdir}/Michael/metrics"

    input:
    path script 
    tuple val(id), path(truth), path(genes), val(pi0), path(guess), val(type)

    output:
    path "metrics_${id}_${pi0}_${type}.txt"

    """
    python3 ${script} ${guess} ${truth} ${genes} ${pi0} ${id} ${type}
    """
}

//Metrics with the null version
process METRICS_NULL{

    container 'metrics:latest'
    publishDir "${params.outdir}/Michael/metrics"

    input:
    path script 
    tuple val(id), path(truth), path(genes), val(pi0), path(guess), val(type)

    output:
    path "metrics_${id}_${pi0}_${type}.txt"

    """
    python3 ${script} ${guess} ${truth} ${genes} ${pi0} ${id} ${type}
    """
}


