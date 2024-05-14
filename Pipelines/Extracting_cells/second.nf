//Here we are going to try and create a file which can read in multiple files from the BEELINE file
//and return something we can test with


//Scripts used
research_con  = file( params.research_con )
add_gene_name = file( params.add_gene_name )
prior_eb      = file( params.prior_eb )
metrics       = file( params.metrics )

//Workflow to import
include { PRIOR_CHANGE as PRIOR_CHANGE_TEST } from './Modules/Workflow'

//Floating Variables
prior_type = "Full"
keep = 0.9
multi_str = "0.01,0.1,1.0,10.0,100.0"
p_val = 0.75
num_genes = 19
data_type = "Curated"
method = "prior_formula_testing"
num_cells = "Not_sure"

//Workflow to test my theory out
workflow {
    //The raw data I am working with
    def groundtruth = Channel.fromPath( "./Data/Research_data/inputs/Curated/GSD/GSD-2000-?/refNetwork.csv" )
    def expressions = Channel.fromPath( "./Data/Research_data/inputs/Curated/GSD/GSD-2000-?/ExpressionData.csv" )
    def multi = Channel.of(0.01,0.1,1,10,100)
    //Credit to this:
    //https://nextflow-io.github.io/patterns/create-key-to-combine-channels/
    expressions
        .map { [it.toString().split('2000-')[1].split('/ExpressionData.csv')[0], it] }
        .set { expression_key }
    groundtruth
        .map { [it.toString().split('2000-')[1].split('/refNetwork.csv')[0], it] }
        .set { groundtruth_key }
    ch_inputs = expression_key
        .join(groundtruth_key, by:0)
        .combine(multi)
    //Convert into something useful
    CONVERT(
        research_con,
        ch_inputs
    )

    //Run Empirical Bayes
    PRIOR_CHANGE_TEST(
        CONVERT.out.tuple_converter,
        add_gene_name,
        prior_type,
        prior_eb,
        keep,
        p_val,
    )

    //Now use these files to get metrics
    ch_output = METRICS(
        metrics,
        PRIOR_CHANGE_TEST.out
    )
}


//Converts what we need into two files which can be now used
process CONVERT {

    publishDir "${params.outdir}/Converted"

    input:
    path script 
    tuple val(id), path(data), path(truth), val(multi)

    output:
    tuple val(id), path("expression${id}_${multi}.txt"), path("ground_truth${id}_${multi}.csv"), 
        val(multi) emit: tuple_converter

    """
    python3 ${script} ${data} ${truth} ${id} ${multi}
    """
}

process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics"

    input:
    path script
    tuple val(id), path(original), path(given)

    output:
    tuple val(Matthew), val(AUROC)

    script:
    Matthew = int
    AUROC = int

    """
    #!/usr/bin/env python3

    import numpy as np
    import sklearn.metrics as metrics
    import sys
    

    #Firstly import the data needed
    given = np.genfromtxt(results/$id/$given, delimiter = ",")
    original = np.genfromtxt(results/$id/$original, delimiter = ",")
    original = [x for row in original for x in row]
    given = [x for row in given for x in row]

    if max(given) == 2:
        $Matthew = None
        $AUROC = None
    else:
        $Matthew = metrics.matthews_corrcoef(original,given)
        fpr, tpr, threshold = metrics.roc_curve(original,given)
        $AUROC = metrics.auc(fpr,tpr)

    """
}


//Process which echos a bunch of files, a way to test channel inputs
// process ECHO {
//     input:
//     tuple val(id), path(file1), path(file2)

//     """
//     echo ${file1} ${file2}
//     """
// }

//We now convert this into a workkable EB thing
// process EB_CONVERT{

//     container 'nlnet_convert:latest'

//     input:
//     path script
//     path infile

//     output:
//     path "formatted_data.txt" emit: data

//     script:

//     """
//     python3 ${script} ${infile}
//     """    

// }