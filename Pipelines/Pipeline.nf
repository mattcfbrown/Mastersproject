#!/usr/bin/env nextflow

params.count = "$projectDir/Data/Test2.txt"
params.outdir = "results"


workflow {
    test_ch = NLNET(params.count)
    test_ch.view()
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


process NLNET {

    container 'docker://brownmc/nlnet-image:latest'

    input:
    path count
    path script

    output:
    path 'matrix.txt'

    script:
    """
    Rscript ${script}

    """
}
