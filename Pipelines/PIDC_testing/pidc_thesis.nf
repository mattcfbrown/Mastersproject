//This file here's only goal is to run the same methods as before for PIDC

//Scripts needed
normalised     = file( params.normalised )
genie_script   = file( params.genie_script )
genie_con      = file( params.genie_con )
nlnet_score    = file( params.nlnet_score )
input_fix      = file( params.input_fix )
pidc_eb        = file( params.pidc_eb )
thesis_met     = file( params.thesis_met )

workflow{
    //Step 1: Read in the data
    //Contains (ID, input, groundtruth, zeroed)

    //The ground truth network
    // def groundtruth_no = Channel.fromPath( "./Data/Original_test.csv" )
    def groundtruth_no = Channel.fromPath( "../../Directed_groundtruth/truth_beeline.csv" )
    //The inputs
    // def input_expressions_no = Channel.fromPath( "../Empirical_Bayes_testing/Data/25genes_*.txt" )
    def input_expressions_no = Channel.fromPath( "../GENCI_nextflow_pipeline/Data/BEELINE/Files/data_*.txt" )
    def zero_prior = Channel.fromPath("../../Directed_groundtruth/zeroed_beeline.csv")
    // def zero_prior = Channel.fromPath("./Data/zeroed.csv")
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    input_expressions_no
        .map { [it.toString().split('data_')[1].split('.txt')[0], it] }
        .set { input_key_no }
    no_prior_inputs = input_key_no
        .combine(groundtruth_no)
        .combine(zero_prior)
    
    //Step 2: Normalise the data
    //This results in:
    //(ID, Normalised, groundtruth, zeroed)
    // normal_ch = NORMALISE(
    //     normalised,
    //     no_prior_inputs
    // )

        //Step 3: RUN GENIE3
    //Returns
    //(ID, normalised, ground truth, zeroed, GENIE output)
    genie_ch = GENIE(
        genie_script,
        no_prior_inputs
    )

    //Step 4: Convert GENIE3 output into something more readable
    //Returns
    //(ID, normalised, ground truth, zeroed, normalised GENIE output)
    genie_con_ch = GENIE_CONVERSION(
        genie_con,
        genie_ch
    )

    //Step 5: Run NLNET
    //Returns
    //(ID, normalised, ground truth, zeroed, normalised GENIE output, nlnet scores)
    // nlnet_ch = NLNET_SCORES(
    //     nlnet_score,
    //     genie_con_ch
    // )

    //Step 6: Run PIDC EB
    convert_ch = EB_CONVERT(
        input_fix,
        genie_con_ch
    )

    output_ch = EB_RUN(
        pidc_eb,
        0.75,
        convert_ch
    )

    METRICS(
        thesis_met,
        output_ch
    )

}

//Normalise the data
process NORMALISE {

    container 'metrics:latest'
    publishDir "${params.outdir}/thesis/normalised"

    input:
    path script 
    tuple val(id), path(input), path(truth), val(zero)

    output:
    tuple val(id), path("normalised_file_${id}.txt"), path(truth), val(zero)

    script:

    """
    python3 ${script} ${input} ${id}
    """   
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
    tuple val(id), path(input), path(truth), path(zero)

    output:
    tuple val(id), path(input), path(truth), path(zero), path("genie${id}.txt")

    script:

    """
    Rscript ${genie_script} ${input} ${id} > genie${id}.txt
    """

}

//Here we will deal with formatting of the output
//This will normalise the data into something which can be used for the output
process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/thesis/genie", mode: 'copy'

    input:
    path script
    tuple val(id), path(input), path(truth), path(zero), path(genie_output)

    output:
    tuple val(id), path(input), path(truth), path(zero), path("matrix_genie_normalised_${id}.csv")

    script:
    """
    python3 ${script} ${genie_output} ${id}
    """
}

//Returns NLENET scores
process NLNET_SCORES {

    container 'nlnet_score:latest'
    publishDir "${params.outdir}/nlnet", mode: 'copy'

    input:
    path script
    tuple val(id), path(input), path(truth), path(zero), path(genie_output)

    output:
    tuple val(id), path(input), path(truth), path(zero), path(genie_output), path("nlnet_scores.csv")

    script:
    """
    Rscript ${script} ${input}
    """

}

//Convert to EB readable format
process EB_CONVERT{

    publishDir "${params.outdir}/prior_genie/EB_formatted"

    input:
    path script
    tuple val(id), path(input), path(truth), path(zero), path(genie_output)

    output:
    tuple val(id), path("formatted_data_${id}.txt"), path(truth), 
        path(zero), path(genie_output), emit: tuple_converter

    script:

    """
    python3 ${script} ${input} ${id}
    """    

}

//Runs Empirical Bayes
process EB_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/thesis/EB_output/prior"

    input:
    path script
    val p_val
    tuple val(id), path(input), path(truth), val(zero), path(genie_output)

    output:
    tuple val(id), path(truth), path(genie_output),
        path("Eb_fullmatrix_${id}.csv"),
        path("Eb_zeromatrix_${id}.csv"),
        path("Eb_geniematrix_${id}.csv"),
        path("PIDC_${id}.csv")

    script:

    """
    julia ${script} ${input} ${truth} ${zero} ${genie_output} ${id} ${p_val}
    """
}

//Metrics time
process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/thesis/metrics/prior"

    input:
    path script 
    tuple val(id), path(truth), path(genie_output), path(eb_full), path(eb_zero), path(eb_genie), path(pidc)

    output:
    path "metrics_${id}.txt"

    """
    python3 ${script} ${eb_zero} ${eb_full} ${eb_genie} ${genie_output} ${truth} ${id} ${pidc}
    """
}