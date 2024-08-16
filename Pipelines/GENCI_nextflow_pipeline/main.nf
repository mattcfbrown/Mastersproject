//This file will run GENECI and produce outputs
//At the moment I want to make sure it works for one run before I start to increase the complexity of the program

//Files used

geneci_convert = file( params.geneci_convert )
input_fix      = file( params.input_fix )
eb_geneci      = file( params.eb_geneci )
metrics        = file( params.metrics )

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
    //(ID, outputs of the grn, optimised ensemble network, EB input, ground truth)

    eb_input_ch = EB_CONVERT(
        input_fix,
        ch_inputs
    )

    //Step 3: runs EB
    //This will return:
    //(ID, outputs of the grn, optimised ensemble network, ground truth, EB full output, EB top 10 output)
    eb_output = EB_RUN(
        eb_geneci,
        0.9,
        eb_input_ch
    )

    //Now it is time to get some metrics
    METRICS(
        metrics,
        eb_output
    )
}

//This converts the data into something we can use
process GENECI_FORMAT {

    publishDir "${params.outdir}/converted_data/GENECI"
    container 'nlnet_convert:latest'

    input:
    path script
    tuple val(id), path(original_input), path(groundtruth)

    output:
    tuple val(id), path(original_input), path("genci_input_${id}.csv"), path(groundtruth)

    script:
    """
    python3 ${script} ${original_input} ${id}
    """
}

//This gives outputs based on GENECI
process INFERENCE {

    publishDir "${params.outdir}"

    input:
    tuple val(id), path(original_input), path(geneci_input), path(groundtruth)

    output:
    tuple val(id), path(original_input), path(groundtruth), path("inference/${id}_folder")

    script:
    """
    source /Users/mbrown/Desktop/Research/Mastersproject/geneci_python/env/bin/activate

    geneci infer-network --expression-data ${geneci_input} \
                         --technique PIDC --technique CLR --technique GENIE3_RF\
                         --output-dir inference/${id}_folder
    """

}

// This runs the optimisation algorithm needed
process OPTIMISE{

    publishDir "${params.outdir}"

    input:
    tuple val(id), path(original_input), path(groundtruth), path(inference_path)

    output:
    tuple val(id), path(original_input), path(groundtruth), path(inference_path), path("optimised/${id}_opti")

    script:
    """
    source /Users/mbrown/Desktop/Research/Mastersproject/geneci_python/env/bin/activate

    geneci optimize-ensemble --confidence-list ${inference_path}/genci_input_${id}/lists/GRN_PIDC.csv \
                             --confidence-list ${inference_path}/genci_input_${id}/lists/GRN_CLR.csv \
                             --confidence-list ${inference_path}/genci_input_${id}/lists/GRN_GENIE3_RF.csv \
                             --crossover-probability 0.9 --mutation-probability 0.05 --population-size 100 \
                             --num-parents 3 --mutation-strength 0.1 \
                             --num-evaluations 1000 --cut-off-criteria PercLinksWithBestConf --cut-off-value 0.4 \
                             --function Quality \
                             --algorithm GA \
                             --output-dir optimised/${id}_opti

    """
}

//Here we convert the input file into something more EB friendly
process EB_CONVERT{

    publishDir "${params.outdir}/EB/formatted"

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

//This program here runs Empirical Bayes
process EB_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/EB/EB_output"

    input:
    path script
    val p_val
    tuple val(id), path(output_path), path(optimise_path), path(input), path(groundtruth)

    output:
    tuple val(id), path(output_path), path(optimise_path), path(groundtruth), 
        path("EB_${id}_full.csv"), path("EB_${id}_ten.csv")

    script:

    """
    julia ${script} ${input} \
          ${optimise_path}/all/final_list.csv \
          ${optimise_path}/top_10/final_list.csv \
          ${p_val} ${id}
    """
}

//Metrics process
process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics"

    input:
    path script
    tuple val(id), path(output_path), path(optimise_path), path(groundtruth), 
        path(eb_full), path(eb_ten)

    output:
    path "Metrics_${id}.txt"

    script:
    """
    python3 ${script} ${output_path}/lists \
            ${eb_full} ${eb_ten} \
            ${groundtruth} ${id} \
            ${optimise_path}/all/final_list.csv \
            ${optimise_path}/top_10/final_list.csv 
    """
}