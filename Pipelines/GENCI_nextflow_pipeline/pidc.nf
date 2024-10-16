//This file here will explore how the PIDC method works for EB when given some ensemble data

//It will have multiple steps

//Files used
input_fix      = file( params.input_fix )
pidc           = file( params.pidc )
metrics        = file( params.metrics )
pidc_boot      = file( params.pidc_boot )

workflow{

    //Step 1: Read in all required data into a channel
    //Get the optimised directory
    def optimised = channel.of(file("./Data/optimised_no_pidc/").listFiles())
    //Get the output directory
    def grn = channel.of(file("./Data/outputs/").listFiles())
    //Get the ground truth networks
    def groundtruth = Channel.fromPath( "./Data/ground_truth/*.csv" )
    //Get the inputs
    def inputs = Channel.fromPath("./Data/Files/*.txt")
    
    // //Credit to this:
    // //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    optimised
        .map { [it.toString().split('genci_input_')[1].split('_folder')[0], it] }
        .set { optimised_key }
    grn
        .map { [it.toString().split('grn_input_')[1].split('_folder')[0], it] }
        .set { grn_key }
    inputs
        .map { [it.toString().split('data_')[1].split('.txt')[0], it] }
        .set { input_key }
    groundtruth
        .map { [it.toString().split('truth_')[1].split('.csv')[0], it] }
        .set { truth_key }
    ch_inputs = grn_key
        .combine(optimised_key, by: 0)
        .combine(input_key, by: 0)
        .combine(truth_key, by: 0)

    //Step 2: EB input correction
    //This will return:
    //(ID, outputs of the grn, optimised ensemble network, EB input, ground truth)

    eb_input_ch = EB_CONVERT(
        input_fix,
        ch_inputs
    )

    //Step 3: EB running
    //Two parts, top 10 and all data
    //Full
    eb_full_ch = EB_RUN_FULL(
        pidc,
        0.9,
        'full',
        eb_input_ch
    )

    //Top 10
    eb_output_ch = EB_RUN_TEN(
        pidc,
        0.9,
        'ten',
        eb_full_ch
    )
    //Step 4: Metrics
    METRICS(
        metrics,
        eb_output_ch
    )
}

//Processes

//Here we convert the input file into something more EB friendly
process EB_CONVERT{

    publishDir "${params.outdir}/PIDC/EB/formatted"

    //ID, outputs of the grn, optimised ensemble network, input file, ground truth

    input:
    path script
    tuple val(id), path(output_path), path(optimise_path), path(input), path(groundtruth)

    output:
    tuple val(id), path(output_path), path(optimise_path), path("formatted_data.txt"), path(groundtruth)

    script:

    """
    python3 ${script} ${input} > formatted_data.txt
    """    

}

//Full prior information
process EB_RUN_FULL {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/PIDC/EB/EB_output"

    input:
    path script
    val p_val
    val type
    tuple val(id), path(output_path), path(optimise_path), path(input), path(groundtruth)

    output:
    tuple val(id), path(output_path), path(optimise_path), path(input), path(groundtruth), 
        path("EB_${id}_${type}.csv")

    script:

    """
    julia ${script} ${input} \
          ${optimise_path}/all/final_list.csv \
          ${p_val} ${id} ${type}
    """
}

//Top 10 prior information
process EB_RUN_TEN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/PIDC/EB/EB_output"

    input:
    path script
    val p_val
    val type
    tuple val(id), path(output_path), path(optimise_path), path(input), path(groundtruth), path(eb_full)

    output:
    tuple val(id), path(output_path), path(optimise_path), path(groundtruth), 
        path(eb_full), path("EB_${id}_${type}.csv"), path(input)

    script:

    """
    julia ${script} ${input} \
          ${optimise_path}/top_10/final_list.csv \
          ${p_val} ${id} ${type}
    """
}

//Bootstrapping process
process BOOTSTRAPPING{

    container 'eb_bootstrap:latest'
    publishDir "${params.outdir}/PIDC/EB/bootstrapping"

    input:
    path script
    val p_val
    tuple val(id), path(output_path), path(optimise_path), path(groundtruth), 
        path(eb_full), path(eb_ten), path(input)
    
    output:
    tuple val(id), path(output_path), path(optimise_path), 
        path("EB_${id}_bootstrapping.csv"), path("EB_${id}_bootstrapping_zero.csv"),
            path(groundtruth), path(eb_full), path(eb_ten)

    script:
    """
    julia -t 4 ${script} ${input} \
          ${optimise_path}/all/final_list.csv \
          ${p_val} \
          ${id} \
    """
}

//Metrics process
process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/PIDC/metrics"

    input:
    path script
    tuple val(id), path(output_path), path(optimise_path), 
        path(eb_boot), path(eb_boot_zero),
            path(groundtruth), path(eb_full), path(eb_ten)

    output:
    path "Metrics_${id}.txt"

    script:
    """
    python3 ${script} ${output_path}/lists \
            ${eb_full} ${eb_ten} \
            ${groundtruth} ${id} \
            ${optimise_path}/all/final_list.csv \
            ${optimise_path}/top_10/final_list.csv \
            ${eb_boot} \
            ${eb_boot_zero} ${eb_zero} \
    """
}

