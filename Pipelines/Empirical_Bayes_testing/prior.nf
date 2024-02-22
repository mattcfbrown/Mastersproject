//This will run multiple tests

//All the scripts needed
file_correct = file( params.input_fix )
bayes_script = file( params.Emp_Bayes )
metric_EB    = file( params.metrics_EB )
genie3_eb    = file( params.genie3_eb )
genie3_con   = file( params.genie3_con )

//All the datasets needed
test_data    = file( params.test_data )
gene_read25  = file( params.gene_read25 )
gene_orig25  = file( params.gene_orig25 )
zero_info    = file( params.no_prior )

//All floating variables
p_value   = 0.9
num_genes = 25

//Now get the packages we wish to work with
include { EMPIRICAL_BAYES as NO_PRIORS } from './modules'
include { EMPIRICAL_BAYES as FULL_PRIORS } from './modules'
include { EMPIRICAL_BAYES as GENIE3_PRIORS } from './modules'
include { GENIE3 as GENIE3 } from './modules'

workflow {

    //Run with no prior information
    no_prior = NO_PRIORS(
        gene_read25,
        zero_info,
        p_value,
        file_correct,
        bayes_script,
        'no_priors'
    )

    //Run with complete prior information
    full_prior = FULL_PRIORS(
        gene_read25,
        gene_orig25,
        p_value,
        file_correct,
        bayes_script,
        'full_priors'
    )

    //Gets genie3 prior information
    genie3_info = GENIE3(
        gene_read25,
        genie3_eb,
        genie3_con,
        num_genes
    )

    //Now performs an empirical Bayes analysis
    genie3_prior = GENIE3_PRIORS(
        gene_read25,
        genie3_info,
        p_value,
        file_correct,
        bayes_script,
        'genie3_priors'
    )

    //Performs a metric analysis
    METRICS(
        metric_EB,
        no_prior,
        full_prior,
        genie3_prior,
        gene_orig25
    )
}

//Here I will perform some metrics on it, using a similar version of the script I developed for my original test
process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics"

    input:
    path script
    path no_prior
    path full_prior
    path genie_prior
    path original

    output:
    path 'ROCplot.pdf'
    path 'matthew.txt'

    script:

    """
    python3 ${script} ${no_prior} ${full_prior} ${genie_prior} ${original}
    """    

}
