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
zero_info10  = file( params.zero_file10 )
gene_read10  = file( params.gene_read10 )
gene_orig10  = file( params.gene_orig10 )

//All floating variables
p_value   = [0.9,0.9,0.9]
num_genes = [10,25,50]
keep = [0.9,1.0,1.0]
type = ["no_priors", "full_priors", "genie3_priors"]

include { PRIOR_WORKFLOW as GENE10 } from './workflow'
include { PRIOR_WORKFLOW as GENE25 } from './workflow'
include { PRIOR_WORKFLOW as GENE50 } from './workflow'

workflow {
    GENE10(
        gene_read10,
        zero_info10,
        gene_orig10,
        p_value[0],
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes[0],
        type,
        keep[0]
    )

    GENE25(
        gene_read25,
        zero_info25,
        gene_orig25,
        p_value[1],
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes[1],
        type,
        keep[1]
    )

    GENE50(
        gene_read50,
        zero_info50,
        gene_orig50,
        p_value[2],
        file_correct,
        bayes_script,
        genie3_eb,
        genie3_con,
        metric_EB,
        num_genes[2],
        type,
        keep[2]
    )
}