//This nextflow script aims to test and see how cell numbers impact the results
//I will use the 25 gene samples
//This will require an update to the metrics function

//All the scripts needed
file_correct = file( params.input_fix )
bayes_script = file( params.Emp_Bayes )
metric_EB    = file( params.metrics_EB )
genie3_eb    = file( params.genie3_eb )
genie3_con   = file( params.genie3_con )

//All the datasets needed
gene_read25_cells100   = file( params.g_r25c100 )            //100 Cells
gene_read25_cells250   = file( params.g_r25c250 )            //250 Cells
gene_read25_cells500   = file( params.g_r25c500 )            //500 Cells
gene_read25_cells1000  = file( params.g_r25c1000 )           //1000 Cells
gene_read25_cells2000  = file( params.g_r25c2000 )           //2000 Cells
zero_info25            = file( params.no_prior )             //This is the zero file
cells_test_orig        = file( params.gene_orig25 )          //This is the original file 

//Import all the workflows we need
include { PRIOR_WORKFLOW as CELL100 } from './workflow'
include { PRIOR_WORKFLOW as CELL250 } from './workflow'
include { PRIOR_WORKFLOW as CELL500 } from './workflow'
include { PRIOR_WORKFLOW as CELL1000 } from './workflow'
include { PRIOR_WORKFLOW as CELL2000 } from './workflow'

//All the Important floating variables
p_value   = 0.9
num_genes = 25
keep = [1.0,1.0,0.9,0.9,0.9]
type = ["no_priors", "full_priors", "genie3_priors"]
num_cells = [100,250,500,1000,2000]
dist = ["gam","gam","normal","normal","normal"]
null_type = ["MI", "MI", "PUC", "PUC", "PUC"]


//Workflow
workflow {

    CELL100(
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
        null_type[0]
    )

    CELL250(
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
        null_type[1]
    )
    
    CELL500(
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
        null_type[2]
    )

    CELL1000(
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
        null_type[3]
    )

    CELL2000(
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
        num_cells[4],
        dist[4],
        null_type[4]
    )        
}