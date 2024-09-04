//This will essentially be the same workflow as before, except with the alternate prior information
//I just need to write a function for EB

//This file will run GENECI and produce outputs
//At the moment I want to make sure it works for one run before I start to increase the complexity of the program

//Files used

geneci_convert = file( params.geneci_convert )
input_fix      = file( params.input_fix )
eb_geneci      = file( params.eb_genci_prior )
metrics        = file( params.metrics )
bootstrap      = file( params.prior_boot )

workflow{

    //Step 1: Load in all the data which I will use
    //This returns a Tuple with the following:
    //(ID, outputs of the grn, optimised ensemble network, input file, ground truth)

    //Get the optimised directory
    def optimised = channel.of(file("./Data/optimised/").listFiles())
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
    //Step 2: The data needs to be put into the acceptable EB format
    //This will return:
    //(ID, outputs of the grn, optimised ensemble network, EB input, original_input, ground truth)

    eb_input_ch = EB_CONVERT(
        input_fix,
        ch_inputs
    )

    //Step 3: Time to Run some bootstrapping methods
    //This will return:
    //(ID, outputs of the grn, optimised ensemble network, EB input, Eb boot, ground truth)
    boot_output_ch = BOOTSTRAPPING(
        bootstrap,
        0.9,
        eb_input_ch
    )


    //Step 4: runs EB
    //This will return:
    //(ID, outputs of the grn, optimised ensemble network, ground truth, EB boot,EB full output, EB top 10 output)
    eb_output = EB_RUN(
        eb_geneci,
        0.9,
        boot_output_ch
    )



    //Step 5: Now it is time to get some metrics
    METRICS(
        metrics,
        eb_output
    )
}


//Here we convert the input file into something more EB friendly
process EB_CONVERT{

    publishDir "${params.outdir}/prior/EB/formatted"

    //ID, outputs of the grn, optimised ensemble network, input file, ground truth

    input:
    path script
    tuple val(id), path(output_path), path(optimise_path), path(input), path(groundtruth)

    output:
    tuple val(id), path(output_path), path(optimise_path), path("formatted_data.txt"), 
        path(input), path(groundtruth)

    script:

    """
    python3 ${script} ${input} > formatted_data.txt
    """    

}

//This program here runs Empirical Bayes
process EB_RUN {

    container 'eb_bootstrap:latest'
    publishDir "${params.outdir}/prior/EB/EB_output"

    input:
    path script
    val p_val
        tuple val(id), path(output_path), path(optimise_path), path(input), 
        path(eb_boot), path(eb_boot_zero), path(groundtruth)

    output:
    tuple val(id), path(output_path), path(optimise_path), path(groundtruth), path(eb_boot), path(eb_boot_zero),
        path("EB_${id}_full.csv"), path("EB_${id}_ten.csv"), path("EB_${id}_five.csv"), path("EB_${id}_zero.csv")

    script:

    """
    julia ${script} ${input} \
          ${optimise_path}/all/final_list.csv \
          ${optimise_path}/top_10/final_list.csv \
          ${optimise_path}/top_5/final_list.csv \
          ${p_val} ${id}
    """
}

//Metrics process
process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/prior/metrics"

    input:
    path script
    tuple val(id), path(output_path), path(optimise_path), path(groundtruth), path(eb_boot), path(eb_boot_zero),
        path(eb_full), path(eb_ten), path(eb_five), path(eb_zero)

    output:
    path "Metrics_${id}.txt"

    script:
    """
    python3 ${script} ${output_path}/lists \
            ${eb_full} ${eb_ten} \
            ${groundtruth} ${id} \
            ${optimise_path}/all/final_list.csv \
            ${optimise_path}/top_10/final_list.csv \
            ${eb_boot} ${eb_five} \
            ${optimise_path}/top_5/final_list.csv \
            ${eb_boot_zero} ${eb_zero} \
    """
}

//Bootstrapping process
process BOOTSTRAPPING{

    container 'eb_bootstrap:latest'
    publishDir "${params.outdir}/prior/EB/bootstrapping"

    input:
    path script
    val p_val
    tuple val(id), path(output_path), path(optimise_path), path(eb_input), 
        path(input), path(groundtruth)
    
    output:
    tuple val(id), path(output_path), path(optimise_path), path(eb_input), 
        path("EB_${id}_bootstrapping.csv"), path("EB_${id}_bootstrapping_zero.csv"),
            path(groundtruth)

    script:
    """
    julia -t 4 ${script} ${input} \
          ${optimise_path}/all/final_list.csv \
          ${p_val} \
          ${id} \
    """
}