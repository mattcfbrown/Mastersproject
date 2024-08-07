//This file will run GENECI and produce outputs
//At the moment I want to make sure it works for one run before I start to increase the complexity of the program

//Files used

geneci_convert = file( params.geneci_convert )

workflow{

    //Step 1: Load in all the data which I will use
    //This returns a Tuple with the following:
    //(ID, input, Ground truth)

    //Ground truth networks
    def groundtruth = Channel.fromPath( "../Ensemble_testing/Data/Original_test.csv" )
    // def groundtruth = Channel.fromPath( "./Data/truth1.csv")
    //The training data
    def inputs = Channel.fromPath( "../Ensemble_testing/Data/25genes_1000cells.txt" )
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


    //Step 2: The data needs to be put into the acceptable GENECI format
    //This will return:
    //(ID, original input, converted input, ground truth)

    format_ch = GENECI_FORMAT(
        geneci_convert,
        ch_inputs
    )

    //Step 3: Run Inference methods
    //This output format will be determined later

    inference_ch = INFERENCE(
        format_ch
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
    container 'geneci_run:latest'

    input:
    tuple val(id), path(original_input), path(geneci_input), path(groundtruth)

    output:
    tuple val(id), path(original_input), path(groundtruth)

    script:
    """
    geneci infer-network --expression-data ${geneci_input} --technique PIDC --output-dir ${params.outdir}
    """

}