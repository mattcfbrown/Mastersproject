//This here will take in several different methods as ensemble priors

//Files used
pearson        = file( params.pearson )
input_fix      = file( params.input_fix )
info_measure   = file( params.info_measure )
eb_methods     = file( params.eb_methods )
metrics_meth   = file( params.metrics_meth )

workflow{
    //The general make up will look as follows

    //Step 1: Read in all the data which is needed
    //This returns a Tuple with the following:
    //(ID, input, Ground truth)

    //Ground truth networks
    def groundtruth = Channel.fromPath( "./Data/Original_test.csv" )
    // def groundtruth = Channel.fromPath( "./Data/truth1.csv")
    //The training data
    def inputs = Channel.fromPath( "./Data/25genes_*cells.txt" )
    // def inputs = Channel.fromPath( "./Data/training*.txt" )
    //The testing data
    //The values of pi0 we will be testing for 2.2 = 0.9, 2.944 = 0.95
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    inputs
        .map { [it.toString().split('25genes_')[1].split('cells.txt')[0], it] }
        // .map { [it.toString().split('training')[1].split('txt')[0], it] }
        .set { expressiontrain_key }
    ch_inputs = expressiontrain_key
        .combine(groundtruth)

    //Step 2: Perform a Pearson and Spearman correlation (plus normalised)
    //This returns the tuple:
    //(ID, input, truth, pearson weight, pearson guess, spearman weight, spearman guess)

    pearson_ch = CORRELATION(
        pearson,
        ch_inputs
    )

    //Step 3: Format the input file ready to be used in the Empirical Bayes programs
    //This returns the tuple:
    //(ID, EB_input, truth, pearson weight, pearson guess, spearman weight, spearman guess)
    converted_ch = EB_CONVERT(
        input_fix,
        pearson_ch
    )

    //Step 4: Run Mutual Information (plus normalised)
    //This returns the tuple:
    //(ID, EB_input, truth, pearson weight, pearson guess, spearman weight, spearman guess, MI weight, MI guess)
    mutual_ch = INFORMATION_MEASURES(
        info_measure,
        converted_ch
    )

    //Step 5: Run Empirical Bayes 
    //This returns the tuple:
    //(ID, EB_input, truth, pearson eb, pearson guess, spearman eb, spearman guess, MI eb, MI guess, zero eb)

    eb_ch = EB_RUN(
        eb_methods,
        0.9,
        mutual_ch
    )

    //Step 6: Get the metrcis needed
    METRICS(
        metrics_meth,
        eb_ch
    )

}

//Methods which I will use:
//Pearson and spearman correlation
process CORRELATION {

    container 'pearson:latest'
    publishDir "${params.outdir}/methods/Pearson"

    input:
    path script
    tuple val(id), path(input), path(truth)

    output:
    tuple val(id), path(input), path(truth), path("pearson_weight_${id}.csv"), path("pearson_guess_${id}.csv"),
        path("spear_weight_${id}.csv"), path("spear_guess_${id}.csv")


    """
    Rscript ${script} ${input} ${id}
    """
}

//Convert file readable by Information measures and EB
process EB_CONVERT{

    publishDir "${params.outdir}/methods/Information_measures/formatted"

    input:
    path script
    tuple val(id), path(input), path(truth), path(pearson_prior), path(pearson_guess), path(spear_weight),
        path(spear_guess)

    output:
    tuple val(id), path("formatted_data.txt"), path(truth), path(pearson_prior), path(pearson_guess),
        path(spear_weight), path(spear_guess)

    script:

    """
    python3 ${script} ${input} > formatted_data.txt
    """    

}

//Run Mutual Information
process INFORMATION_MEASURES {

    container "empiricalbayes:latest"
    publishDir "${params.outdir}/methods/Information_measures/outputs"

    input:
    path script
    tuple val(id), path(input), path(truth), path(pearson_prior), path(pearson_guess),
        path(spear_weight), path(spear_guess)

    output:
    tuple val(id), path(input), path(truth), path(pearson_prior), path(pearson_guess),
        path(spear_weight), path(spear_guess),
            path("MI_scores_${id}.csv"), path("MI_matrix_${id}.csv"), emit: tuple_converter

    script:

    """
    julia ${script} ${input} ${id}
    """   
}

//Run the Empirical Bayes framework
process EB_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/methods/EB_output"

    input:
    path script
    val p_val
    tuple val(id), path(input), path(truth), path(pearson_prior), path(pearson_guess),
        path(spear_weight), path(spear_guess),
            path(mi_weight), path(mi_guess)

    output:
    tuple val(id), path(input), path(truth), path("Eb_pearson_${id}.csv"), path(pearson_guess),
        path("Eb_spear_${id}.csv"), path(spear_guess),
            path("Eb_MI_${id}.csv"), path(mi_guess), path("Eb_zero_${id}.csv")

    script:

    """
    julia ${script} ${input} ${pearson_prior} ${spear_weight} ${mi_weight} ${p_val} ${id}
    """
}

//Metrics
process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/methods/metrics"

    input:
    path script
    tuple val(id), path(input), path(truth), path(eb_pearson), path(pearson_guess),
        path(eb_spear), path(spear_guess),
            path(eb_mi), path(mi_guess), path(eb_zero)

    output:
    path "matthew_${id}.txt"

    script:
    """
    python3 ${script} ${truth} ${eb_pearson} ${pearson_guess} ${eb_spear} ${spear_guess} ${eb_mi} ${mi_guess} ${eb_zero} ${id}
    """
}
