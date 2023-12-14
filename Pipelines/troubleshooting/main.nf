

// declare nlnet rscript IS a file (so we can pass it correctly to process)
nlnet_rscript  = file( params.nlnet_rscript )
NI_script      = file( params.NI_script )
EB_script      = file( params.EB_script )
count_matrix   = file( params.count_data )
SCENIC_test    = file( params.test_SC)
SCENIC_matrix  = file( params.matrix)
Trans_fact     = file ( params.TFs)

workflow {
    //PIDC METHOD FOR NETWORK INFERENCE
    //INFORMATION_MEASURES(
    //    count_matrix,
    //    NI_script,
    //    0.15
    //).view()
    //Empirical Bayes method for network inference
    //EMPIRICAL_BAYES(
    //    count_matrix,
    //    EB_script
    //).view()
    //nlnet method for network inference
    //NLNET(
    //    count_matrix,
    //    nlnet_rscript
    //).view()
    SCENIC(
        SCENIC_matrix,
        Trans_fact
    ).view()    

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

    container 'camaralab/scenic:latest'
    
    input:
    path grn
    path TF
    
    output:
    path "outfile.loom", emit: outfile
    
    script:

    """
    
    pyscenic grn ${grn} ${TF} > outfile.loom

    """
}