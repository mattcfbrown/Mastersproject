//This here is a program which is able to test how effective prior usage is for different softwares


//All the scripts needed
file_correct = file( params.input_fix )

//All the datasets needed
test_data    = file( params.test_data )
gene_read25  = file( params.gene_read25 )

//Workflow
workflow {

    EB_CONVERT(
        gene_read25,
        file_correct
    ).view()

}

//This process converts the reads into something empirical Bayes can use
process EB_CONVERT{

    publishDir "${params.outdir}/Empirical_Bayes"

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

