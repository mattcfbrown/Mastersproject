

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


workflow {
    //PIDC METHOD FOR NETWORK INFERENCE
    //INFORMATION_MEASURES(
    //    Original,
    //    NI_script,
    //    0.15,
    //    Network
    //).view()
    //Empirical Bayes method for network inference
    //EMPIRICAL_BAYES(
    //    count_EB,
    //    EB_script
    //).view()

    //nlnet method for network inference
    //Firstly get the output needed
    output_nlnet = NLNET(
        count_matrix,
        nlnet_rscript
    )
    //Now get it into a readable version
    nlnet_matrix = NLNET_CONVERSION (
        output_nlnet,
        nlenet_con_scr,
        25
    )
    //Now compare to the actual matrix
    //NOTE, I WILL CONTINUE TO UPDATE THIS COMMAND AS I GET ALL THE RESULTS IN
    METRICS(
        nlnet_matrix,
        original_mat,
        metric_program
    ).view()

    //SCENIC(
    //    SCENIC_matrix,
    //    Trans_fact
    //).view()    

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

process INFORMATION_MEASURES {

    container 'networkinference:latest'
    publishDir "${params.outdir}/Information_Measures"

    input:
    path infile
    path jlscript
    val threshold
    path network

    output:
    path 'outfile_NI.txt'

    script:

    """

    julia ${jlscript} ${infile} ${threshold} ${network} > outfile_NI.txt

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
    path original_matrix
    path metric_script

    output:
    path 'ROCplot.pdf'

    script:

    """
    python3 ${metric_script} ${nlnet_matrix} ${original_matrix} > ROCplot.pdf
    """    

}