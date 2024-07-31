//This program is a catch all for attempting to derive some ensemble based methods

//Paths used
nlnet_script   = file( params.nlnet )
genie_met_con  = file( params.genie_met_con )
genie_script   = file( params.genie_script )
input_fix      = file( params.input_fix )
info_measure   = file( params.info_measure )
ensemble       = file( params.ensemble )
eb_run         = file( params.eb_run )
metrics        = file( params.metrics )
normalisation  = file( params.normalisation )
normalise_eb   = file( params.normalise_eb )
ensemble_mat   = file( params.ensemble_mat )

//NOTE: Statistics includes: Normalised confidence scores AND a plot of what these scores will look like
    //This also includes a prediction about how the network should look like

//Rough plan

workflow{
    //Step 1: Put all needed data into tuples (This will be for the first one initially, add more once the pipeline works)

    //This returns a Tuple with the following:
    //(ID, Testing data, Ground truth)

    //Ground truth networks
    // def groundtruth = Channel.fromPath( "./Data/Original_test.csv" )
    def groundtruth = Channel.fromPath( "./Data/truth1.csv")
    //The training data
    // def inputs = Channel.fromPath( "./Data/25genes_*cells.txt" )
    def inputs = Channel.fromPath( "../Actual_data_tests/Research_data/inputs/Curated/HSC/HSC-2000-?-50/ExpressionData.csv" )
    //The testing data
    //The values of pi0 we will be testing for 2.2 = 0.9, 2.944 = 0.95
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    inputs
        // .map { [it.toString().split('25genes_')[1].split('cells.txt')[0], it] }
        .map { [it.toString().split('2000-')[1].split('-50/ExpressionData.csv')[0], it] }
        .set { expressiontrain_key }
    ch_inputs = expressiontrain_key
        .combine(groundtruth)


    //Step 2: Run NLNET, get some NLNET statistics

    //Shut down for the moment
    // nlnet_output = NLNET(
    //     nlnet_script,
    //     ch_inputs
    // )


    //Step 3: Run Genie, get some GENIE statistics
    //This will create a tuple with the following values
    // (ID, converted testing data, truth matrix, prior matrix, genie prediction matrix)

    //Run Genie
    genie3_ch = GENIE(
        genie_script,
        ch_inputs
    )

    //Convert it all
    prior_ch = GENIE_CONVERSION(
        genie_met_con,
        genie3_ch
    )

    //Step 4: Run Julia, using 3/4 of the inference methods, get statistics
    //This returns the tuple:
    //(ID, input, truth, genie prior, genie guess, mi prior, mi guess, puc prior, puc guess, clr prior,
    // clr guess)

    //Convert into a format EB can work with
    prior_ch = EB_CONVERT(
        input_fix,
        prior_ch
    )

    //We now want to output some information measures
    ensemble_ch = INFORMATION_MEASURES(
        info_measure,
        prior_ch
    )


    //Step 5: Ensemble all the data (Weighted averages)
    ensemble_prior = ENSEMBLE(
        ensemble,
        ensemble_ch
    )

    //Step 6: Try to see if the ensemble method has a good threshold
    ensemble_prior_2 = ENSEMBLE_MATRIX(
        ensemble_mat,
        ensemble_prior
    )

    //Step 7: Run Empirical Bayes
    //This returns the following Tuple:
    //(ID, ground truth, EB_guess, Genie_guess, MI_guess, PUC_guess, CLR_guess)
    eb_output_ch = EB_RUN(
        eb_run,
        0.8,
        ensemble_prior_2
    )


    //Step 8: Metrics
    //Compare it's performance to all other methods, we will see what gets spit out in the end
    METRICS(
        metrics,
        eb_output_ch
    )

    //I also need to make some graphs:
    NORMALISATION(
        normalisation,
        ensemble_ch
    )

    EB_NORMALISATION(
        normalise_eb,
        eb_output_ch
    )    
}

//Processes to be run
//NLNET
process NLNET {

    container 'nlnet_score:latest'
    publishDir "${params.outdir}/nlnet", mode: 'copy'

    input:
    path nlnet_script
    tuple val(id), path(input), path(truth)

    output:
    tuple val(id), path(input), path(truth), path("nlnet_scores_${id}.csv"), path("nlnet_${id}.csv")

    script:
    """
    Rscript ${nlnet_script} ${input} ${id}
    """

}

//Here we run GENIE3 on the testing data
//We do no normalising here
process GENIE{

    cpus 4
    memory '4 GB'
    container 'genie3:latest'
    publishDir "${params.outdir}/genie/genie_outputs"

    input:
    path genie_script
    tuple val(id), path(input), path(truth)

    output:
    tuple val(id), path(input), path(truth), path("genie${id}.txt")

    script:

    """
    Rscript ${genie_script} ${input} ${id} > genie${id}.txt
    """

}

//Here we will deal with formatting of the output
//We also need to have an outcome for the matrix, therefore a threshold value shall be used
process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/genie/priors", mode: 'copy'

    input:
    path script
    tuple val(id), path(input), path(truth), path(genie_output)

    output:
    tuple val(id), path(input), path(truth), path("matrix_genie_eb_${id}.csv"), path("genie_predict_${id}.csv")

    script:
    """
    python3 ${script} ${genie_output} ${id}
    """
}

//Here we convert the file into something useable
process EB_CONVERT{

    publishDir "${params.outdir}/Information_measures/formatted"

    input:
    path script
    tuple val(id), path(input), path(truth), path(prior), path(genie_guess)

    output:
    tuple val(id), path('formatted_data.txt'), path(truth), 
        path(prior), path(genie_guess), emit: tuple_converter

    script:

    """
    python3 ${script} ${input} > formatted_data.txt
    """    

}

process INFORMATION_MEASURES {

    container "empiricalbayes:latest"
    publishDir "${params.outdir}/Information_measures/outputs"

    input:
    path script
    tuple val(id), path(input), path(truth), path(prior), path(genie_guess)

    output:
    tuple val(id), path(input), path(truth), path(prior), path(genie_guess),
        path("MI_scores_${id}.csv"), path("MI_matrix_${id}.csv"), 
            path("PUC_scores_${id}.csv"), path("PUC_matrix_${id}.csv"),
                path("CLR_scores_${id}.csv"), path("CLR_matrix_${id}.csv"), emit: tuple_converter

    script:

    """
    julia ${script} ${input} ${id}
    """   
}

process ENSEMBLE {

    container 'metrics:latest'
    publishDir "${params.outdir}/ensemble"

    input:
    path script
    tuple val(id), path(input), path(truth), path(genie_prior), path(genie_guess),
        path(mi_score), path(mi_guess), 
            path(puc_score), path(puc_guess),
                path(clr_score), path(clr_guess)
    
    output:
    tuple val(id), path(input), path(truth),
        path(genie_guess), path(mi_guess), path(puc_guess), path(clr_guess),
            path("ensemble_${id}.csv")

    script:
    """
    python3 ${script} ${id} ${truth} ${genie_prior} ${mi_guess} ${puc_guess} ${clr_guess}
    """
}

process EB_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/EB_output"

    input:
    path script
    val p_val
    tuple val(id), path(input), path(truth),
        path(genie_guess), path(mi_guess), path(puc_guess), path(clr_guess),
            path(ensemble_prior), path(ensemble_guess)

    output:
    tuple val(id), path(truth), path("Eb_matrix_${id}.csv"),
        path(genie_guess), path(mi_guess), path(puc_guess), path(clr_guess),
            path(ensemble_prior), path(ensemble_guess)

    script:

    """
    julia ${script} ${input} ${ensemble_prior} ${p_val} ${id}
    """
}

process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics"

    input:
    path script
    tuple val(id), path(truth), path(eb_guess),
        path(genie_guess), path(mi_guess), path(puc_guess), path(clr_guess),
            path(ensemble_prior), path(ensemble_guess)

    output:
    path "matthew_${id}.txt"

    script:
    """
    python3 ${script} ${truth} ${eb_guess} ${genie_guess} ${mi_guess} ${puc_guess} ${clr_guess} ${ensemble_guess} ${id}
    """
}

process NORMALISATION {
    
    container 'metrics:latest'
    publishDir "${params.outdir}/metrics/normalisation_graphs"

    input:
    path script
    tuple val(id), path(input), path(truth), path(genie_score), path(genie_guess),
        path(mi_score), path(mi_guess), 
            path(puc_score), path(puc_guess),
                path(clr_score), path(clr_guess)

    output:
    path "normalisation_${id}_Genie.pdf"
    path "normalisation_${id}_MI.pdf"
    path "normalisation_${id}_PUC.pdf"
    path "normalisation_${id}_CLR.pdf"

    script:
    """
    python3 ${script} ${truth} ${id} ${genie_guess} ${mi_guess} ${puc_guess} ${clr_guess} ${genie_score} ${mi_score} ${puc_score} ${clr_score}
    """
}

process EB_NORMALISATION {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics/normalisation_graphs"   
    
    input:
    path script
    tuple val(id), path(truth), path(eb_guess),
        path(genie_guess), path(mi_guess), path(puc_guess), path(clr_guess),
            path(ensemble_prior), path(ensemble_guess)

    output:
    path "normalisation_${id}_EB.pdf"

    script:
    """
    python3 ${script} ${truth} ${id} ${ensemble_prior} ${eb_guess}
    """
}

process ENSEMBLE_MATRIX {

    container 'metrics:latest'
    publishDir "${params.outdir}/ensemble"
    
    input:
    path script
    tuple val(id), path(input), path(truth),
        path(genie_guess), path(mi_guess), path(puc_guess), path(clr_guess),
            path(ensemble_prior)
    
    output:
    tuple val(id), path(input), path(truth),
        path(genie_guess), path(mi_guess), path(puc_guess), path(clr_guess),
            path(ensemble_prior), path("ensemble_${id}_matrix.csv")
    
    script:
    """
    python3 ${script} ${ensemble_prior} ${id}
    """
}