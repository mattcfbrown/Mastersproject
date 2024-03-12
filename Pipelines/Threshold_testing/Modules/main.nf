//This here will be where we can put all the workflows in one place



//Getting the packages
//This is the mutual information one
include { INFORMATION_MEASURES as MI_THRESHOLD1 } from './Workflow'
include { INFORMATION_MEASURES as MI_THRESHOLD2 } from './Workflow'
include { INFORMATION_MEASURES as MI_THRESHOLD3 } from './Workflow'
include { INFORMATION_MEASURES as MI_THRESHOLD4 } from './Workflow'
include { INFORMATION_MEASURES as MI_THRESHOLD5 } from './Workflow'
//This is the GENIE ones
include { GENIE3 as GENIE3_THRESHOLD1 } from './Workflow'
include { GENIE3 as GENIE3_THRESHOLD2 } from './Workflow'
include { GENIE3 as GENIE3_THRESHOLD3 } from './Workflow'
include { GENIE3 as GENIE3_THRESHOLD4 } from './Workflow'
include { GENIE3 as GENIE3_THRESHOLD5 } from './Workflow'


//Workflows
//Information measures
workflow THRESH_INFORMATION_MEASURES{
    take:
    current                     //Our original count matrix
    file_correct               //Puts our count matrix into a workable form
    NI_script                 //The script which runs Information measures
    NI_con_scr               //This script converts the Information measures output into something which can be measured
    threshold               //This is the threshold used
    num_genes              //The number of genes in the dataset
    original              //Original script
    type                 //The type of data used
    method              //The method employed
    num_cells          //Number of cells in the data
    metric            //Metrics script needed
    threshold_str    //A string with the threshold values


    main:
    output1 = MI_THRESHOLD1(
        current,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold[0],
        num_genes
    )

    output2 = MI_THRESHOLD2(
        current,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold[1],
        num_genes
    )

    output3 = MI_THRESHOLD3(
        current,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold[2],
        num_genes
    )

    output4 = MI_THRESHOLD4(
        current,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold[3],
        num_genes
    )

    output5 = MI_THRESHOLD5(
        current,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold[4],
        num_genes
    )

    METRIC_THRESHOLD(
        metric,
        output1,
        output2,
        output3,
        output4,
        output5,
        original,
        threshold_str,
        type,
        method,
        num_cells,
        num_genes
    )
}
//GENIE3
workflow THRESH_GENIE3{
    take:
    current                    //Our original count matrix
    genie_script              //The script which runs GENIE3
    genie_con                //This script converts the GENIE3 output into something which can be measured
    threshold               //This is the threshold used
    num_genes              //The number of genes in the dataset
    original              //Original script
    type                 //The type of data used
    method              //The method employed
    num_cells          //Number of cells in the data
    metric            //Metrics script needed
    threshold_str    //A string with the threshold values

    main:
    output1 = GENIE3_THRESHOLD1(
        current,
        genie_script,
        genie_con,
        threshold[0],
        num_genes
    )

    output2 = GENIE3_THRESHOLD2(
        current,
        genie_script,
        genie_con,
        threshold[1],
        num_genes
    )

    output3 = GENIE3_THRESHOLD3(
        current,
        genie_script,
        genie_con,
        threshold[2],
        num_genes
    )

    output4 = GENIE3_THRESHOLD4(
        current,
        genie_script,
        genie_con,
        threshold[3],
        num_genes
    )

    output5 = GENIE3_THRESHOLD5(
        current,
        genie_script,
        genie_con,
        threshold[4],
        num_genes
    )

    METRIC_THRESHOLD(
        metric,
        output1,
        output2,
        output3,
        output4,
        output5,
        original,
        threshold_str,
        type,
        method,
        num_cells,
        num_genes
    )
}


//Metrics process
process METRIC_THRESHOLD {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics"

    input:
    path metric
    path threshold1
    path threshold2
    path threshold3    
    path threshold4
    path threshold5
    path original
    val threshold_str
    val type
    val method
    val num_cells
    val num_genes



    output:
    path "ROCplot_${num_genes}_${num_cells}_${type}_${method}.pdf"
    path "matthew_${num_genes}_${num_cells}_${type}_${method}.txt"

    script:

    """
    python3 ${metric} ${threshold1} ${threshold2} ${threshold3} ${threshold4} ${threshold5} ${original} ${threshold_str} ${type} ${method} ${num_cells}
    """   
}