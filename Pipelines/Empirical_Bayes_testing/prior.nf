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
num_genes = 25
type = ["no_priors", "full_priors", "genie3_priors"]

include { GENE25 as PRIOR_WORKFLOW } from './workflow'

workflow {
    GENE25(
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
}