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
zero_info25  = file( params.no_prior )
zero_info50  = file( params.zero_file50 )
gene_read50  = file( params.gene_read50 )
gene_orig50  = file( params.gene_orig50 )

//All floating variables
p_value   = 0.9
num_genes = 50
type = ["no_priors", "full_priors", "genie3_priors"]

//Now get the packages we wish to work with
include { EMPIRICAL_BAYES as NO_PRIORS } from './modules'
include { EMPIRICAL_BAYES as FULL_PRIORS } from './modules'
include { EMPIRICAL_BAYES as GENIE3_PRIORS } from './modules'
include { GENIE3 as GENIE3 } from './modules'

workflow PRIOR_WORKFLOW {
    take:
    reads                                      //The RNA sequencing reads
    null_prior                                //CSV file full of zeros, no info provided
    known_network                            //The known GRN
    p_value                                 //The p value we test against
    eb_input_fix                           //This fixes the input given to the empirical Bayes function
    bayes_script                          //Runs Empirical Bayes
    eb_genie                             //Runs GENIE3 in the context of Empirical Bayes
    genie_con                           //Runs a conversion script which converts the output of GENIE3 into something workable
    metric_eb                          //Empirical Bayes metrics script
    num_genes                         //Number of genes in the RNA seq file
    type                             //A list which holds all the string data

    main:
    //Run with no prior information
    no_prior = NO_PRIORS(
        reads,
        null_prior,
        p_value,
        eb_input_fix,
        bayes_script,
        type[0]
    )
    //Run with complete prior information
    full_prior = FULL_PRIORS(
        reads,
        known_network,
        p_value,
        eb_input_fix,
        bayes_script,
        type[1]
    )
    //Gets genie3 prior information
    genie3_info = GENIE3(
        reads,
        eb_genie,
        genie3_con,
        num_genes
    )
    //Now performs an empirical Bayes analysis
    genie3_prior = GENIE3_PRIORS(
        reads,
        genie3_info,
        p_value,
        eb_input_fix,
        bayes_script,
        type[2]
    )
    //Performs a metric analysis
    METRICS(
        metric_eb,
        no_prior,
        full_prior,
        genie3_prior,
        known_network,
        num_genes
    )
}

workflow {

    PRIOR_WORKFLOW(
        gene_read25,
        zero_info25,
        gene_orig25,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type        
    )

    // //Run with no prior information
    // no_prior = NO_PRIORS(
    //     gene_read50,
    //     zero_info50,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     'no_priors'
    // )

    // //Run with complete prior information
    // full_prior = FULL_PRIORS(
    //     gene_read50,
    //     gene_orig50,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     'full_priors'
    // )

    // //Gets genie3 prior information
    // genie3_info = GENIE3(
    //     gene_read50,
    //     genie3_eb,
    //     genie3_con,
    //     num_genes
    // )

    // //Now performs an empirical Bayes analysis
    // genie3_prior = GENIE3_PRIORS(
    //     gene_read50,
    //     genie3_info,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     'genie3_priors'
    // )

    // //Performs a metric analysis
    // METRICS(
    //     metric_EB,
    //     no_prior,
    //     full_prior,
    //     genie3_prior,
    //     gene_orig50,
    //     num_genes
    // )
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
    val num_genes 

    output:
    path "ROCplot_${num_genes}.pdf"
    path "matthew_${num_genes}.txt"

    script:

    """
    python3 ${script} ${no_prior} ${full_prior} ${genie_prior} ${original} ${num_genes}
    """    

}
