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
beeline_g1  = file( params.beeline_g1 )            
beeline_g2  = file( params.beeline_g2 )       
beeline_g3  = file( params.beeline_g3 )       
beeline_g4  = file( params.beeline_g4 )       
beeline_g5  = file( params.beeline_g5 )                

beeline_m1  = file( params.beeline_m1 )
beeline_m2  = file( params.beeline_m2 )       
beeline_m3  = file( params.beeline_m3 )       
beeline_m4  = file( params.beeline_m4 )       
beeline_m5  = file( params.beeline_m5 )       

zero_info25            = file( params.zeroed_bee )             //This is the zero file
cells_test_orig        = file( params.beeline_orig )          //This is the original file 


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

//All the Important floating variables
p_value   = 0.9990909091
num_genes = 11
keep = [0.8,0.8,0.8,0.8,0.8]
type = ["no_priors", "full_priors", "genie3_priors", "nlnet_priors"]
num_cells = ['test1','test2','test3','test4','test5']
dist = ["gam","gam","normal","normal","normal"]
null_type = ["MI", "MI", "PUC","PUC", "PUC"]
threshold = 0.1


//Workflow
workflow {

    CELL100(
        // normalised,
        beeline_g1,
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
        // normalised,
        beeline_m1,
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
        0.8,
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
        // normalised,
        beeline_g2,
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
        // normalised,
        beeline_m2,
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
        // normalised,
        beeline_g3,
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
        // normalised,
        beeline_m3,
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
        // normalised,
        beeline_g4,
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
        // normalised,
        beeline_m4,
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


    CELL2000(
        // normalised,
        beeline_g5,
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
        // normalised,
        beeline_m5,
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
        0.8,
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
}
