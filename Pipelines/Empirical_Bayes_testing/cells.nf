//This nextflow script aims to test and see how cell numbers impact the results
//I will use the 25 gene samples
//This will require an update to the metrics function

//All the scripts needed
file_correct   = file( params.input_fix )
bayes_script   = file( params.Emp_Bayes )
metric_EB      = file( params.metrics_EB )
genie3_eb      = file( params.genie3_eb )
genie3_con     = file( params.genie3_con )
nlnet_score    = file( params.nlnet_score )
genie3_nor     = file( params.genie3_nor )
genie3_nor_cor = file( params.genie3_nor_cor )
NI_script      = file( params.NI_script )
NI_con_scr     = file( params.NI_convert )
normalised     = file( params.normalised)

//All the datasets needed
gene_read25_cells100   = file( params.g_r25c100 )            //100 Cells
gene_read25_cells250   = file( params.g_r25c250 )            //250 Cells
gene_read25_cells500   = file( params.g_r25c500 )            //500 Cells
gene_read25_cells1000  = file( params.g_r25c1000 )           //1000 Cells
gene_read25_cells2000  = file( params.g_r25c2000 )           //2000 Cells

messy100               = file( params.messy100 )             //Messy data 100 cells
messy250               = file( params.messy250 )             //Messy data 250 cells
messy500               = file( params.messy500 )             //Messy data 500 cells
messy1000              = file( params.messy1000 )            //Messy data 1000 cells
messy2000              = file( params.messy2000 )            //Messy data 2000 cells

zero_info25_old        = file( params.no_prior )             //This is the zero file
zero_info25            = file( params.half )                 //This is the 0.5 prior file, it has be written with this name so I do not have to edit the workflows below
cells_test_orig        = file( params.gene_orig25 )          //This is the original file 

clean_test1            = file( params.clean_test1 )
clean_test2            = file( params.clean_test2 )
clean_test3            = file( params.clean_test3 )
clean_test4            = file( params.clean_test4 )
messy_test1            = file( params.messy_test1 )
messy_test2            = file( params.messy_test2 )
messy_test3            = file( params.messy_test3 )
messy_test4            = file( params.messy_test4 )
original_test1         = file( params.original_test1 )
original_test2         = file( params.original_test2 )
original_test3         = file( params.original_test3 )
original_test4         = file( params.original_test4 )

//Import all the workflows we need
//Clean Workflows
include { PRIOR_WORKFLOW as CELL100 } from './workflow'
include { PRIOR_WORKFLOW as CELL250 } from './workflow'
include { PRIOR_WORKFLOW as CELL500 } from './workflow'
include { PRIOR_WORKFLOW as CELL1000 } from './workflow'
include { PRIOR_WORKFLOW as CELL2000 } from './workflow'

//Messy workflows
include { PRIOR_WORKFLOW as MESSY100 } from './workflow'
include { PRIOR_WORKFLOW as MESSY250 } from './workflow'
include { PRIOR_WORKFLOW as MESSY500 } from './workflow'
include { PRIOR_WORKFLOW as MESSY1000 } from './workflow'
include { PRIOR_WORKFLOW as MESSY2000 } from './workflow'


//2000 CELL workflows
include { PRIOR_WORKFLOW as CLEAN_TEST_1 } from './workflow'
include { PRIOR_WORKFLOW as CLEAN_TEST_2 } from './workflow'
include { PRIOR_WORKFLOW as CLEAN_TEST_3 } from './workflow'
include { PRIOR_WORKFLOW as CLEAN_TEST_4 } from './workflow'
include { PRIOR_WORKFLOW as MESSY_TEST_1 } from './workflow'
include { PRIOR_WORKFLOW as MESSY_TEST_2 } from './workflow'
include { PRIOR_WORKFLOW as MESSY_TEST_3 } from './workflow'
include { PRIOR_WORKFLOW as MESSY_TEST_4 } from './workflow'


//All the Important floating variables
p_value   = 0.75
num_genes = 25
keep = [1.0,1.0,0.9,0.9,0.9]
type = ["no_priors", "full_priors", "genie3_priors", "nlnet_priors"]
num_cells = [100,250,500,1000,2000]
dist = ["gam","gam","normal","normal","normal"]
null_type = ["MI", "MI", "PUC","PUC", "PUC"]
threshold = 0.1


//Workflow
workflow {

    CELL100(
        normalised,
        gene_read25_cells100,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[0],
        num_cells[0],
        dist[0],
        null_type[0],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    MESSY100(
        normalised,
        messy100,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        0.7,
        "Messy_100",
        dist[2],
        null_type[2],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    CELL250(
        normalised,
        gene_read25_cells250,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[1],
        num_cells[1],
        dist[1],
        null_type[1],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    MESSY250(
        normalised,
        messy250,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[1],
        "Messy_250",
        dist[1],
        null_type[1],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )
    
    CELL500(
        normalised,
        gene_read25_cells500,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[2],
        num_cells[2],
        dist[2],
        null_type[2],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )
    
    MESSY500(
        normalised,
        messy500,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[3],
        "Messy_500",
        dist[3],
        null_type[3],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    CELL1000(
        normalised,
        gene_read25_cells1000,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[3],
        num_cells[3],
        dist[3],
        null_type[3],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    MESSY1000(
        normalised,
        messy1000,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[3],
        "Messy_1000",
        dist[3],
        null_type[3],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    //This is the 2000 cells test

    CELL2000(
        normalised,
        gene_read25_cells2000,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[4],
        "2000_Test_0",
        dist[4],
        null_type[4],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    MESSY2000(
        normalised,
        messy2000,
        zero_info25,
        cells_test_orig,
        p_value,
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes,
        type,
        keep[4],
        "Messy_2000_Test_0",
        dist[4],
        null_type[4],
        nlnet_score,
        NI_script,
        NI_con_scr,
        threshold,
        genie3_nor,
        genie3_nor_cor
    )

    // CLEAN_TEST_1(
    //     clean_test1,
    //     zero_info25,
    //     original_test1,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "2000_Test_1",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )

    // CLEAN_TEST_2(
    //     clean_test2,
    //     zero_info25,
    //     original_test2,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "2000_Test_2",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )  

    // CLEAN_TEST_3(
    //     clean_test3,
    //     zero_info25,
    //     original_test3,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "2000_Test_3",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )

    // CLEAN_TEST_4(
    //     clean_test4,
    //     zero_info25,
    //     original_test4,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "2000_Test_4",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )

    // MESSY_TEST_1(
    //     messy_test1,
    //     zero_info25,
    //     original_test1,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "Messy_2000_Test_1",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )

    // MESSY_TEST_2(
    //     messy_test2,
    //     zero_info25,
    //     original_test2,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "Messy_2000_Test_2",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )      

    // MESSY_TEST_3(
    //     messy_test3,
    //     zero_info25,
    //     original_test3,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "Messy_2000_Test_3",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )          

    // MESSY_TEST_4(
    //     messy_test4,
    //     zero_info25,
    //     original_test4,
    //     p_value,
    //     file_correct,
    //     bayes_script,
    //     genie3_eb,
    //     genie3_con,
    //     metric_EB,
    //     num_genes,
    //     type,
    //     keep[4],
    //     "Messy_2000_Test_4",
    //     dist[4],
    //     null_type[4],
    //     nlnet_score,
    //     NI_script,
    //     NI_con_scr,
    //     threshold,
    //     genie3_nor,
    //     genie3_nor_cor
    // )                                 
}
