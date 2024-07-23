// This file will be running GENIE3 priors, using the SERGIO data
//The aim here is to get some pretty figures to use 

//The scripts we will use
genie_met_con  = file( params.genie_met_con )
genie_script   = file( params.genie_script )
genie_con_gene = file( params.genie_con_gene )
eb_and_met     = file( params.genie_and_met )
input_fix      = file( params.input_fix )
genie_metric   = file( params.genie_metric )

//A tuple of values which will be used for the Empirical Bayes run
eb_tuple = tuple(0.75, 0.80)

workflow{

    //STEP 1: Get all data needed into the appropriate tuples

    //We will get a tuple with the following structure:
    //(ID, training dataset, testing dataset, ground truth network)

    //Ground truth networks
    def groundtruth = Channel.fromPath( "./Data/Original_test.csv" )
    //The training data
    def training_expressions = Channel.fromPath( "./Data/25genes_*cells.txt" )
    //The testing data
    def testing_expressions = Channel.fromPath( "./Data/25genes_*cells_messy.txt" )
    //The values of pi0 we will be testing for 2.2 = 0.9, 2.944 = 0.95
    def pi0 = Channel.of(2.2,2.944)    
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    training_expressions
        .map { [it.toString().split('25genes_')[1].split('cells.txt')[0], it] }
        .set { expressiontrain_key }
    testing_expressions
        .map { [it.toString().split('25genes_')[1].split('cells_messy.txt')[0], it] }
        .set { expressiontest_key }
    ch_inputs = expressiontrain_key
        .join(expressiontest_key, by:0)
        .combine(groundtruth)

    //STEP 2: Run GENIE3 on the 50 file (using no priors)
    
    //This will create a tuple with the following values
    // (ID, testing data, truth matrix, genie3 outputs)
    genie3_ch = GENIE(
        genie_script,
        ch_inputs
    )

    //STEP 3: Normalise the output and put reference network in a matrix

    //This will create a tuple with the following values
    // (ID, converted testing data, truth matrix, prior matrix, genie prediction matrix)

    prior_ch = GENIE_CONVERSION(
        genie_met_con,
        genie3_ch
    )


    //STEP 4: With the reference network, run EB again, this time on the messy data using clean prior

    //Firstly convert the training data into something worthwhile
    prior_ch = EB_CONVERT(
        input_fix,
        prior_ch
    )


    //Run Empirical Bayes
    //This will return a tuple of the following nature
    // (Id, true network, genie's network, the EB network, data from running EB, gene_lists, prior information)
    eb_ch = EB_RUN(
        eb_and_met,
        eb_tuple,
        prior_ch
    )

    //Step 5: We get the metrics
    METRICS(
        genie_metric,
        eb_ch
    )
}


//Here we run GENIE3 on the testing data
//We do no normalising here
process GENIE{

    cpus 4
    memory '4 GB'
    container 'genie3:latest'
    publishDir "${params.outdir}/geniemetrics/genie_outputs"

    input:
    path genie_script
    tuple val(id), path(training), path(testing), path(truth)

    output:
    tuple val(id), path(testing), path(truth), path("genie${id}.txt")

    script:

    """
    Rscript ${genie_script} ${training} ${id} > genie${id}.txt
    """

}

//Here we will deal with formatting of the output
//We also need to have an outcome for the matrix, therefore a threshold value shall be used
process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/geniemetrics/priors", mode: 'copy'

    input:
    path script
    tuple val(id), path(testing), path(truth), path(genie_output)

    output:
    tuple val(id), path(testing), path(truth), path("matrix_genie_eb_${id}.csv"), path("genie_predict_${id}.csv")

    script:
    """
    python3 ${script} ${genie_output} ${id}
    """
}

//Here we convert the file into something useable
process EB_CONVERT{

    publishDir "${params.outdir}/geniemetrics/EB_formatted"

    input:
    path script
    tuple val(id), path(testing), path(truth), path(prior), path(genie_guess)

    output:
    tuple val(id), path('formatted_data.txt'), path(truth), 
        path(prior), path(genie_guess), emit: tuple_converter

    script:

    """
    python3 ${script} ${testing} > formatted_data.txt
    """    

}

//Here we run Empirical Bayes
process EB_RUN {

    container 'changing_priors:latest'
    publishDir "${params.outdir}/geniemetrics/EB_output"

    input:
    path script
    tuple val(p_val), val(keep)
    tuple val(id), path(testing), path(truth), path(prior), path(genie_guess)

    output:
    tuple val(id), path(truth), path(genie_guess), path("EB_matrix_${id}.csv"), path("data_${id}.csv"),
        path("gene_list_${id}.txt"), path(prior)

    script:

    """
    julia ${script} ${testing} ${prior} ${p_val} ${id}
    """
}

process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/geniemetrics/metrics"

    input:
    path script
    tuple val(id), path(truth), path(genie_guess), path(eb_guess), path(data), path(genes), path(prior)

    output:
    path("metrics_${id}.txt")
    path("Agreed_${id}.pdf")
    path("Genie_${id}.pdf")
    path("Prior_${id}.pdf")
    path("Incorrect_Genie_${id}.pdf")
    path("Incorrect_Prior_${id}.pdf")
    path("Null_Mixture_Agreed_${id}.pdf") 
    path("Null_Mixture_Genie_${id}.pdf") 
    path("Null_Mixture_Prior_${id}.pdf") 
    path("Null_Mixture_Incorrect_Genie_${id}.pdf") 
    path("Null_Mixture_Incorrect_prior_${id}.pdf")   
    path("mixture_agree_${id}.pdf")   
    path("mixture_genie_${id}.pdf")   
    path("mixture_prior_${id}.pdf")   
    path("mixture_incorrect_prior_${id}.pdf")   
    path("mixture_incorrect_genie_${id}.pdf")     
    
    script:

    """
    python3 ${script} ${genie_guess} ${eb_guess} ${truth} ${data} ${genes} ${prior} ${id}

    """

}