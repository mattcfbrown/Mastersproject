

// declare nlnet rscript IS a file (so we can pass it correctly to process)
nlnet_rscript  = file( params.nlnet_rscript )
NI_script      = file( params.NI_script )
EB_script      = file( params.EB_script )
count_matrix   = file( params.count_data )
count_EB       = file( params.EB_count )
SCENIC_test    = file( params.test_SC)
SCENIC_matrix  = file( params.matrix)
Trans_fact     = file( params.TFs)
Network        = file( params.network)
Original       = file( params.original)
nlenet_con_scr = file( params.nlnet_convert)
original_mat   = file( params.original_mat)
metric_program = file( params.metric)
NI_con_scr     = file( params.NI_convert)
file_correct   = file( params.input_fix)
genie_script   = file( params.genie3_rscript)



// workflow SUBWORKFLOW {
//     // the inputs of this subworkflow
//     take: 
//     path count_data
//     path script

//     // yuour nextflow processes
//     run:
//     MYPROCESS(
//         count_data,
//         script
//     )

//     // the outputs of this subworkflow
//     emit: 
//     path MYPROCESS.outfile 
// }



workflow {
    //PIDC METHOD FOR NETWORK INFERENCE
    //STEP 1: Put the data into a readable format
    //corrected_matrix = INPUT_INFORMATION_MEASURES(
    //    count_matrix,
    //    file_correct
    //)

    //STEP 2: Perform the algorithm
    //output_NI = INFORMATION_MEASURES(
    //    corrected_matrix,
    //    NI_script,
    //    0.15
    //)

    //STEP 3: Convert it into a useable form
    //ni_matrix = NI_CONVERSION(
    //    output_NI,
    //    NI_con_scr,
    //    25
    //)

    //Empirical Bayes method for network inference
    //EMPIRICAL_BAYES(
    //    count_EB,
    //    EB_script
    //).view()

    //nlnet method for network inference
    //Firstly get the output needed
    //output_nlnet = NLNET(
    //    count_matrix,
    //    nlnet_rscript
    //)
    //Now get it into a readable version
    //nlnet_matrix = NLNET_CONVERSION (
    //    output_nlnet,
    //    nlenet_con_scr,
    //    25
    //)
    //Now compare to the actual matrix
    //NOTE, I WILL CONTINUE TO UPDATE THIS COMMAND AS I GET ALL THE RESULTS IN
    //METRICS(
    //    nlnet_matrix,
    //    ni_matrix,
    //    original_mat,
    //    metric_program
    //).view()

    //SCENIC(
    //    SCENIC_matrix,
    //    Trans_fact
    //).view()

    //GENIE3 Method
    GENIE3(
        count_matrix,
        genie_script,
        0.1
    ).view()    

}

process NLNET {

    container 'nlnet:latest'
    publishDir "${params.outdir}/nlnet"

    input:
    path infile
    path rscript

    output:
    path 'outfile_nlnet.txt'

    script:
    """
    Rscript ${rscript} ${infile} > outfile_nlnet.txt
    """
}

process NLNET_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/nlnet"

    input:
    path output_nlnet
    path nlnet_converter
    val num_genes

    output:
    path 'matrix_nlnet.csv'

    script:
    """
    python3 ${nlnet_converter} ${output_nlnet} ${num_genes} > matrix_nlnet.csv
    """

}

process INPUT_INFORMATION_MEASURES {

    publishDir "${params.outdir}/Information_Measures"

    input:
    path infile
    path script

    output:
    path 'formatted_data.txt'

    script:

    """
    python3 ${script} ${infile} > formatted_data.txt
    """
}

process INFORMATION_MEASURES {

    container 'networkinference:latest'
    publishDir "${params.outdir}/Information_Measures"

    input:
    path infile
    path jlscript
    val threshold

    output:
    path 'outfile_NI.txt'

    script:

    """

    julia ${jlscript} ${infile} ${threshold} > outfile_NI.txt

    """

}

process NI_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/Information_Measures"

    input:
    path output_NI
    path NI_converter
    val num_genes

    output:
    path 'matrix_NI.csv'

    script:
    """
    python3 ${NI_converter} ${output_NI} ${num_genes} > matrix_NI.csv
    """
}

process GENIE3 {

    container 'genie3:latest'
    publishDir "${params.outdir}/genie3"

    input:
    path infile
    path genie_script
    val threshold

    output:
    path 'genie.txt'

    script:

    """
    Rscript ${genie_script} ${infile} ${threshold} > genie.txt
    """

}

process EMPIRICAL_BAYES {
    
    container 'empiricalbayes:latest'

    input:
    path infile
    path jlscript

    output:
    path 'outfile_EB.txt'

    script:

    """

    julia ${jlscript} ${infile} > outfile_EB.txt

    """
}

process SCENIC {

    container 'aertslab/pyscenic:0.9.18'
    publishDir "${params.outdir}/SCENIC"
    
    input:
    path grn
    path TF
    
    output:
    path "3A_seurat_RNA_pySCENIC_GRN_adjacencies.csv", emit: outfile
    
    script:

    """
    
    pyscenic grn ${grn} ${TF} -o expr_mat.adjacencies.tsv


    """
}

process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics"

    input:
    path nlnet_matrix
    path ni_matrix
    path original_matrix
    path metric_script

    output:
    path 'ROCplot.pdf'

    script:

    """
    python3 ${metric_script} ${nlnet_matrix} ${ni_matrix} ${original_matrix} > ROCplot.pdf
    """    

}