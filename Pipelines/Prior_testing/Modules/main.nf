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
    prior_type                     //This is the type of prior we are giving it
    keep                          //The proportion of values we are keeping when running EB
    w0                           //The parameter used to calculate the prior information
    num_genes                   //The number of genes in the dataset
    original                   //Original script
    type                      //The type of data used
    method                   //The method employed
    num_cells               //Number of cells in the data
    metric                 //Metrics script needed
    threshold_str         //A string with the threshold values
    extra                //extra string so I don't have to rewrite my threshold function


    main:
    threshold_1 = w0_TEST_1(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        w0[0],
        extra
    )

    threshold_2 = w0_TEST_2(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        w0[1],
        extra
    )

    threshold_3 = w0_TEST_3(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        w0[2],
        extra
    )

    threshold_4 = w0_TEST_4(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        w0[3],
        extra
    )

    threshold_5 = w0_TEST_5(
        reads,
        priors,
        p_val,
        conversion_script,
        eb_script,
        type,
        keep,
        w0[4],
        extra
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
    reads                                           //The data we are analysing
    conversion_script                              //Converts the data into something useable
    prior_type                                    //The type of prior information we are using
    prior_script                                 //The script which runs out EB test
    priors                                      //The prior information used
    to_keep                                    //The percentage of the values to keep
    multi                                     //How we are affecting the priors
    num_genes                                //The number of genes in the dataset
    original                                //Original script
    data_type                              //The type of data used
    method                                //The method employed
    num_cells                            //Number of cells in the data
    metric                              //Metrics script needed
    threshold_str                      //A string with the threshold values
    p_val                             //Threshold value used for 


    main:
    PRIOR_CHANGE_TEST_1(
        reads,
        conversion_script,
        prior_type,
        prior_script,
        priors,
        to_keep,
        multi[0],
        p_val
    )

    PRIOR_CHANGE_TEST_2(
        reads,
        conversion_script,
        prior_type,
        prior_script,
        priors,
        to_keep,
        multi[1],
        p_val
    )

    PRIOR_CHANGE_TEST_3(
        reads,
        conversion_script,
        prior_type,
        prior_script,
        priors,
        to_keep,
        multi[2],
        p_val
    )

    PRIOR_CHANGE_TEST_4(
        reads,
        conversion_script,
        prior_type,
        prior_script,
        priors,
        to_keep,
        multi[3],
        p_val
    )

    PRIOR_CHANGE_TEST_5(
        reads,
        conversion_script,
        prior_type,
        prior_script,
        priors,
        to_keep,
        multi[4],
        p_val
    )

    METRIC_THRESHOLD_PRIOR(
        metric,
        PRIOR_CHANGE_TEST_1.out,
        PRIOR_CHANGE_TEST_2.out,
        PRIOR_CHANGE_TEST_3.out,
        PRIOR_CHANGE_TEST_4.out,
        PRIOR_CHANGE_TEST_5.out,
        original,
        threshold_str,
        data_type,
        method,
        num_cells,
        num_genes
    ) 

}


