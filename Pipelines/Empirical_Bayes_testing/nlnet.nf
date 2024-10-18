//The main purpose of this test is to examine where and why the NLNET is failing
//May also work on GENIE3 priors

genie3_eb      = file( params.genie3_eb )
genie3_con     = file( params.genie3_con )
data  = file( params.g_r25c1000 )
file_correct   = file( params.input_fix )
bayes_script   = file( params.Emp_Bayes )

workflow{

    genie = GENIE3_RUN(
        data,
        genie3_eb
    )

    scores = GENIE_CONVERSION(
        genie,
        genie3_con,
        25
    )

    eb_data = EB_CONVERT(
        data,
        file_correct,
        1000
    )

    EMPIRICAL_BAYES_RUN(
        eb_data,
        bayes_script,
        scores,
        0.9998,
        'bug_fixing',
        0.8,
        "normal",
        "PUC"
    )
}


process GENIE3_RUN {

    cpus 4
    memory '4 GB'
    container 'genie3:latest'
    publishDir "${params.outdir}/bug_fixing/genie3"

    input:
    path infile
    path genie_script

    output:
    path 'genie.txt'

    script:

    """
    Rscript ${genie_script} ${infile} > genie.txt
    """

} 

process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/bug_fixing/genie3", mode: 'copy'

    input:
    path output_genie
    path genie_converter
    val num_genes

    output:
    path 'matrix_genie_eb.csv'

    script:
    """
    python3 ${genie_converter} ${output_genie} ${num_genes}
    """
}

//This process converts the reads into something empirical Bayes can use
process EB_CONVERT{

    publishDir "${params.outdir}/bug_fixing/eb_convert"

    input:
    path infile
    path script
    val num_cells

    output:
    path "formatted_data_${num_cells}.txt"

    script:

    """
    python3 ${script} ${infile} ${num_cells}
    """    

}

process EMPIRICAL_BAYES_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/bug_fixing/EB"

    input:
    path data
    path script
    path priors
    val p_val
    val type
    val keep
    val gam_or_norm
    val inference

    output:
    path "Eb_matrix_${type}.csv"

    script:

    """
    julia ${script} ${data} ${p_val} ${priors} ${type} ${keep} ${gam_or_norm} ${inference}
    """
   
}