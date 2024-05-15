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
type = ["no_priors", "full_priors", "genie3_priors", "nlnet_priors"]

//Now get the packages we wish to work with
include { EMPIRICAL_BAYES as NO_PRIORS } from './modules'
include { EMPIRICAL_BAYES as FULL_PRIORS } from './modules'
include { EMPIRICAL_BAYES as GENIE3_PRIORS } from './modules'
include { EMPIRICAL_BAYES as NLNET_PRIORS} from './modules'
include { GENIE3_EB as GENIE3_EB } from './modules'
include { NLNET as NLNET } from './modules'
//These are the workflows for normal processes
include { GENIE3 as GENIE3 } from './modules'
include { INFORMATION_MEASURES as INFORMATION_MEASURES} from './modules'

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
    keep                            //Tells us the proportion to keep
    num_cells                      //The number of cells in the data
    gam_or_norm                   //This states whether to use Gamma or Normal fitting    
    inference                    //This says what type of inference technique to use
    nlnet_script                //This is the nlnet scores script
    NI_script                  //The script which runs Information measures
    NI_con_scr                //This script converts the Information measures output into something which can be measured
    threshold                //This is the threshold used    
    genie_script            //The script which runs GENIE3
    genie_con_nor          //This script converts the GENIE3 output into something which can be measured

    main:
    //Run with no prior information
    no_prior = NO_PRIORS(
        reads,
        null_prior,
        p_value,
        eb_input_fix,
        bayes_script,
        type[0],
        keep,
        gam_or_norm,
        inference
    )
    //Run with complete prior information
    full_prior = FULL_PRIORS(
        reads,
        known_network,
        p_value,
        eb_input_fix,
        bayes_script,
        type[1],
        keep,
        gam_or_norm,
        inference
    )
    //Gets genie3 prior information
    genie3_info = GENIE3_EB(
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
        type[2],
        keep,
        gam_or_norm,
        inference
    )
    //We now do the same thing for nlnet
    nlnet_info = NLNET(
        reads,
        nlnet_script
    )
    //Now gets the NLNET prior
    nlnet_prior = NLNET_PRIORS(
        reads,
        nlnet_info,
        p_value,
        eb_input_fix,
        bayes_script,
        type[3],
        keep,
        gam_or_norm,
        inference
    )

    genie3_normal = GENIE3(
        reads,
        genie_script,
        genie_con_nor,
        threshold,
        num_genes
    )

    ni_normal = INFORMATION_MEASURES(
        reads,
        eb_input_fix,
        NI_script,
        NI_con_scr,
        threshold,
        num_genes   
    )

    //Performs a metric analysis
    METRICS(
        metric_eb,
        no_prior,
        full_prior,
        genie3_prior,
        nlnet_prior,
        ni_normal,
        genie3_normal,
        known_network,
        num_genes,
        num_cells
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
    path nlnet_prior
    path ni_normal
    path genie3_normal
    path original
    val num_genes
    val num_cells

    output:
    path "ROCplot_${num_genes}_${num_cells}.pdf"
    path "matthew_${num_genes}_${num_cells}.txt"

    script:

    """
    python3 ${script} ${no_prior} ${full_prior} ${genie_prior} ${nlnet_prior} ${genie3_normal} ${ni_normal} ${original} ${num_genes} ${num_cells}
    """    

}

//This is where I get the graphical information pertaining to 
process GRAPHS {

    container 'posterior:latest'                  //I need a new container here
    publishDir "${params.outdir}/metrics/graphs"

    input:
    path script
    path data
    val to_keep
    val p_value
    val data_type
    val num_cells

    output:
    path "Null_Mix_${data_type}_${num_cells}.pdf"
    path "Posterior_${data_type}_${num_cells}.pdf"

    script:

    """
    julia ${script} ${data} ${to_keep} ${p_value} ${data_type} ${num_cells}
    """
}
