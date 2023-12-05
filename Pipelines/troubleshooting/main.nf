

// declare nlnet rscript IS a file (so we can pass it correctly to process)
nlnet_rscript  = file( params.nlnet_rscript )
count_matrix   = file( params.count_data )

workflow {
    NLNET(
        count_matrix,
        nlnet_rscript
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

    input:
    path count

    output:
    path 'matrix.txt'

    script:

    """
    """

}

