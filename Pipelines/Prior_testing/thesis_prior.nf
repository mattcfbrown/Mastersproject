//This will be the outputs which will help me with my thesis

//Scripts used
normalised     = file( params.normalised )
genie_script   = file( params.genie_script )
genie_con      = file( params.genie_con )
input_fix      = file( params.input_fix )
eb_prior       = file( params.eb_prior )
thesis_metrics = file( params.thesis_metrics )
thresh_thesis  = file( params.thresh_thesis )


//Main program here
workflow{

    //Step 1: Read in all the data
    //This produces:
    //(ID, Datafile, ground truth, zeroed)

    //The ground truth network
    def groundtruth_no = Channel.fromPath( "./Data/Original_test.csv" )
    // def groundtruth_no = Channel.fromPath("../../Directed_groundtruth/truth_beeline.csv")
    //The inputs
    def input_expressions_no = Channel.fromPath( "./Data/25genes_*.txt" )
    // def input_expressions_no = Channel.fromPath("../GENCI_nextflow_pipeline/Data/BEELINE/Files/data_*.txt")
    // def zero_prior = Channel.fromPath("../../Directed_groundtruth/zeroed_beeline.csv")
    def zero_prior = Channel.fromPath("../PIDC_testing/Data/zeroed.csv")
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    input_expressions_no
        // .map { [it.toString().split('data_')[1].split('.txt')[0], it] }
        .map { [it.toString().split('25genes_')[1].split('.txt')[0], it] }
        .set { input_key_no }
    inputs_ch = input_key_no
        .combine(groundtruth_no)
        .combine(zero_prior)
    //Step 2: Normalise the data
    //Returns
    //(ID, normalised, ground truth, zeroed)
    // normalised_ch = NORMALISE(
    //     normalised,
    //     inputs_ch
    // )

    //Step 3: RUN GENIE3
    //Returns
    //(ID, normalised, ground truth, zeroed, GENIE output)
    genie_ch = GENIE(
        genie_script,
        inputs_ch
    )

    //Step 4: Convert GENIE3 output into something more readable
    //Returns
    //(ID, normalised, ground truth, zeroed, normalised GENIE output)
    genie_con_ch = GENIE_CONVERSION(
        genie_con,
        genie_ch
    )

    //Step 5: Convert into a readable EB format
    //Returns
    //(ID, EB normalised, ground truth, zeroed, normalised GENIE output)
    eb_readable_ch = EB_CONVERT(
        input_fix,
        genie_con_ch
    )

    //Step 6: Add the different values to the channels
    //Returns
    //(ID, EB normalised, ground truth, zeroed, normalised GENIE output, w0)
    // def multi = Channel.of(1,10,100)
    // def multi = Channel.of(5)
    def pi0 = Channel.of(2.2,2.944)
    pi0_ch = eb_readable_ch
        .combine(pi0)
        // .combine(multi)
    //Step 7: Work with EB
    eb_ch = EB_RUN(
        eb_prior,
        // thresh_thesis,
        0.75,
        pi0_ch
    )

    //Step 8: Metrics
    METRICS(
        thesis_metrics,
        eb_ch
    )

}

//Normalise the data
process NORMALISE {

    container 'metrics:latest'
    publishDir "${params.outdir}/normalised"

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

//Here we convert the file into something useable
process EB_CONVERT{

    publishDir "${params.outdir}/thesis/EB_formatted"

    input:
    path script
    tuple val(id), path(input), path(truth), path(zero), path(genie)

    output:
    tuple val(id), path("formatted_data.txt"), path(truth), 
        path(zero), path(genie), emit: tuple_converter

    script:

    """
    python3 ${script} ${input} ${id} > formatted_data.txt
    """    

}

process EB_RUN {

    container 'eb_bootstrap:latest'
    publishDir "${params.outdir}/thesis/EB_output"

    input:
    path script
    val pval
    tuple val(id), path(input), path(truth), path(zero), path(genie), val(w0)

    output:
    tuple val(id), path(truth), path(genie), val(w0), 
        path("Eb_fullmatrix_${id}_${w0}.csv"),
        path("Eb_geniematrix_${id}_${w0}.csv"),
        path("Eb_zeromatrix_${id}_${w0}.csv")

    script:

    """
    julia ${script} ${input} ${w0} ${truth} ${genie} ${zero} ${pval} ${id}
    """
}

//Metrics time
process METRICS{

    container 'metrics:latest'
    publishDir "${params.outdir}/thesis/metrics/SERGIO"

    input:
    path script 
    tuple val(id), path(truth), path(genie), val(w0), 
        path(ebfull), path(ebgenie), path(ebzero)
    output:
    path "matthew_${id}_${w0}.txt"

    """
    python3 ${script} ${ebzero} ${ebfull} ${ebgenie} ${genie} ${truth} ${id} ${w0}
    """
}