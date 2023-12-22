

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

workflow {
    //PIDC METHOD FOR NETWORK INFERENCE
    //INFORMATION_MEASURES(
    //    count_matrix,
    //    NI_script,
    //    0.15,
    //    Network
    //).view()
    //Empirical Bayes method for network inference
    EMPIRICAL_BAYES(
        count_EB,
        EB_script
    ).view()
    //nlnet method for network inference
    //NLNET(
    //    count_matrix,
    //    nlnet_rscript
    //).view()
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
    path 'Rplots.pdf'

    script:
    """
    Rscript ${rscript} ${infile}
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
    
    input:
    path grn
    path TF
    
    output:
    path "outfile.loom", emit: outfile
    
    script:

    """
    
    pyscenic grn ${grn} ${TF} -o outfile.loom

    """
}