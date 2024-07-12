//This file here is designed to test out how the PIDC method fairs

//Scripts uses
genie_script   = file( params.genie_script )
genie_con_gene = file( params.genie_con_gene )
input_fix      = file( params.input_fix )
running_eb     = file( params.pidc_eb )

//P_value we will use
p_val = 0.9

workflow{
    //Step 1: Get all the files needed into tuples:
    //I need three tuples:
        //Tuple 1: No prior information
        //Tuple 2: Genie3
        //Tuple 3: Full priors
    
    //Tuple 1: No prior information

    //Contains (ID, input, groundtruth, no_prior_identifier, prior information)

    //The ground truth network
    def groundtruth_no = Channel.fromPath( "./Data/Original_test.csv" )
    //The inputs
    def input_expressions_no = Channel.fromPath( "./Data/25genes_*cells.txt" )
    def prior_id_no = Channel.of("No_prior")
    def zero_prior = Channel.fromPath("./Data/zeroed.csv")
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    input_expressions_no
        .map { [it.toString().split('25genes_')[1].split('cells.txt')[0], it] }
        .set { input_key_no }
    no_prior_inputs = input_key_no
        .combine(groundtruth_no)
        .combine(prior_id_no)
        .combine(zero_prior)
    //Tuple 2: GENIE

    //Contains (ID, input, groundtruth, genie_identifier)

    //The ground truth network
    def groundtruth_genie = Channel.fromPath( "./Data/Original_test.csv" )
    //The inputs
    def input_expressions_genie = Channel.fromPath( "./Data/25genes_*cells.txt" )
    def prior_id_genie = Channel.of("genie_prior")   
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    input_expressions_genie
        .map { [it.toString().split('25genes_')[1].split('cells.txt')[0], it] }
        .set { input_key_genie }
    genie_prior_inputs = input_key_genie
        .combine(groundtruth_genie)
        .combine(prior_id_genie)

    //Tuple 3: Full Prior

    //Contains (ID, input, groundtruth, full_identifier, full prior information)

    //The ground truth network
    def groundtruth_full = Channel.fromPath( "./Data/Original_test.csv" )
    //The inputs
    def input_expressions_full = Channel.fromPath( "./Data/25genes_*cells.txt" )
    def prior_id_full = Channel.of("full_prior")
    def full_prior = Channel.fromPath( "./Data/Full_priors.csv" )
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    input_expressions_full
        .map { [it.toString().split('25genes_')[1].split('cells.txt')[0], it] }
        .set { input_key_full }
    full_prior_inputs = input_key_full
        .combine(groundtruth_full)
        .combine(prior_id_full)
        .combine(full_prior)


    //STEP 2: Run GENIE3 on the second tuple group

    //Returns (ID, input ,ground truth, Genie identifier, genie output)

    genie_ch = GENIE(
        genie_script,
        genie_prior_inputs
    )

    //Step 3: Clean up the GENIE3 output
    //Returns (ID, input ,ground truth, identifier, genie prior)

    genie_ch = GENIE_CONVERSION(
        genie_con_gene,
        genie_ch
    )
    
    //Step 4: Combine all three channels
    //Returns (ID, input, ground truth, identifier, prior information)
    combined_channels = no_prior_inputs.concat(genie_ch,full_prior_inputs)


    //Step 5: Run Empirical Bayes
    //This will leave us with (ID, ground truth, prior type, EB output)

    //Firstly convert it into something worthwhile
    combined_channels = EB_CONVERT(
        input_fix,
        combined_channels
    )

    //Secondly run EB
    eb_outputs = EB_RUN(
        running_eb,
        p_val,
        combined_channels
    )

    //Step 6: We get some metrics about it's performance
    

}

//Processes used

process GENIE{

    cpus 4
    memory '4 GB'
    container 'genie3:latest'
    publishDir "${params.outdir}/geniemetrics/genie_outputs"

    input:
    path genie_script
    tuple val(id), path(input), path(truth), val(identifier)

    output:
    tuple val(id), path(input) ,path(truth), val(identifier), path("genie${id}.txt")

    script:

    """
    Rscript ${genie_script} ${input} ${id} > genie${id}.txt
    """

}

process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/prior_genie/priors", mode: 'copy'

    input:
    path script
    tuple val(id), path(input), path(truth), val(identifier), path(genie_output)

    output:
    tuple val(id), path(input), path(truth), val(identifier), path("matrix_genie_normalised${id}.csv")

    script:
    """
    python3 ${script} ${genie_output} ${id}
    """
}

process EB_CONVERT{

    publishDir "${params.outdir}/prior_genie/EB_formatted"

    input:
    path script
    tuple val(id), path(input), path(truth), val(identifier), path(prior)

    output:
    tuple val(id), path("formatted_data_${id}.txt"), path(truth), 
        val(identifier), path(prior), emit: tuple_converter

    script:

    """
    python3 ${script} ${input} ${id}
    """    

}

//Here we run Empirical Bayes
process EB_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/prior_genie/EB_output"

    input:
    path script
    val p_val
    tuple val(id), path(input), path(truth), val(identifier), path(prior)

    output:
    tuple val(id), path(truth), val(identifier), path("Eb_matrix_${id}_${identifier}.csv")

    script:

    """
    julia ${script} ${input} ${prior} ${p_val} ${id} ${identifier}
    """
}

