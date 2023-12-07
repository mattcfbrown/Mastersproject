

// declare nlnet rscript IS a file (so we can pass it correctly to process)
nlnet_rscript  = file( params.nlnet_rscript )
NI_script      = file( params.NI_script )
EB_script      = file( params.EB_script )
count_matrix   = file( params.count_data )

workflow {
    //PIDC METHOD FOR NETWORK INFERENCE
    //INFORMATION_MEASURES(
    //    count_matrix,
    //    NI_script,
    //    0.15
    //).view()
    //Empirical Bayes method for network inference
    EMPIRICAL_BAYES(
        count_matrix,
        EB_script
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

