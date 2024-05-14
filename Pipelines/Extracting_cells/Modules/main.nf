//This is where we put the workflows which are needed


//Getting the packages
//w0 changes
include { EMPIRICAL_BAYES as w0_TEST_1 } from './Workflow'
include { EMPIRICAL_BAYES as w0_TEST_2 } from './Workflow'
include { EMPIRICAL_BAYES as w0_TEST_3 } from './Workflow'
include { EMPIRICAL_BAYES as w0_TEST_4 } from './Workflow'
include { EMPIRICAL_BAYES as w0_TEST_5 } from './Workflow'
//Prior changes
include { PRIOR_CHANGE as PRIOR_CHANGE_TEST_1} from './Workflow'
include { PRIOR_CHANGE as PRIOR_CHANGE_TEST_2} from './Workflow'
include { PRIOR_CHANGE as PRIOR_CHANGE_TEST_3} from './Workflow'
include { PRIOR_CHANGE as PRIOR_CHANGE_TEST_4} from './Workflow'
include { PRIOR_CHANGE as PRIOR_CHANGE_TEST_5} from './Workflow'

include { THRESHOLD as METRIC_THRESHOLD_W0 } from './Workflow'
include { THRESHOLD as METRIC_THRESHOLD_PRIOR } from './Workflow'

//Workflows
//w0 Test
workflow w0_TESTS {
    take:
    reads                               //This is the gene sequencing reads
    priors                             //This is the prior information 
    p_val                             //This is the p value
    conversion_script                //Script which converts the files into something workable
    eb_script                       //This script runs Empirical Bayes
    priot_type                     //This is the type of prior we are giving it
    keep                          //The proportion of values we are keeping when running EB
    gam_or_norm                  //This states whether to use Gamma or Normal fitting
    inference                   //What type of inference to use
    w0                         //The parameter used to calculate the prior information
    num_genes                 //The number of genes in the dataset
    original                 //Original script
    type                    //The type of data used
    method                 //The method employed
    num_cells             //Number of cells in the data
    metric               //Metrics script needed
    threshold_str       //A string with the threshold values


    main:
    threshold_1 = w0_TEST_1(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        gam_or_norm,
        inference,
        w0[0],
    )

    threshold_2 = w0_TEST_2(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        gam_or_norm,
        inference,
        w0[1],
    )

    threshold_3 = w0_TEST_3(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        gam_or_norm,
        inference,
        w0[2],
    )

    threshold_4 = w0_TEST_4(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        gam_or_norm,
        inference,
        w0[3],
    )

    threshold_5 = w0_TEST_5(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        gam_or_norm,
        inference,
        w0[4],
    )

    METRIC_THRESHOLD_W0(
        metric,
        threshold_1,
        threshold_2,
        threshold_3,
        threshold_4,
        threshold_5,
        original,
        threshold_str,
        type,
        method,
        num_cells,
        num_genes
    )
}

//workflow 
workflow PRIOR_CHANGE_TEST {
    take:
    ch_data                                         //Tuple containing: ID, Reads, Prior information
    conversion_script                              //Converts the data into something useable
    prior_type                                    //The type of prior information we are using
    prior_script                                 //The script which runs out EB test
    to_keep                                     //The percentage of the values to keep
    p_val                                      //Threshold value used for 


    main:
    PRIOR_CHANGE_TEST_1(
        ch_data,
        conversion_script,
        prior_type,
        prior_script,
        to_keep,
        p_val
    )
}


