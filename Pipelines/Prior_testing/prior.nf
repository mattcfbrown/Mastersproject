//This workflow aims to see how editing certain properties of the Empirical Bayes prior distribution 
//changes the values we get

//There will be two different tests conducted:
//1: Altering w0
//2: Altering the formula

//Scripts needed
input_fix     = file( params.input_fix )
eb_script     = file( params.Emp_bay )
metrics       = file( params.metrics )
prior_eb      = file( params.prior_eb )
genie_script  = file( params.genie_script )
genie_con     = file( params.genie_con )

//Data needed
reads    = file( params.testing_data )
messy    = file( params.messy_data )
original = file( params.original )

//Floating values
p_val = 0.75
data_type = ["Clean","Messy"]
prior_type = "Full"
method = ["w0_testing","prior_formula_testing"]
keep = 0.9
num_genes = 25
num_cells = 1000
extra = 'extra'

w0 = [-2.0,-1.0,0.0,1.0,2.0]
threshold_str = "-2.0,-1.0,0.0,1.0,2.0"

multi = [0.01,0.1,1.0,10.0,100.0]
multi_str = "0.01,0.1,1.0,10.0,100.0"

//Loading in the workflows needed
include { w0_TESTS as w0_TESTS } from './Modules'
include { PRIOR_CHANGE_TEST as PRIOR_CHANGE_TEST } from './Modules'

//GENIE3 workflow used
workflow GENIE3 {
    //Inputs
    take:
    current                    //Our original count matrix
    genie_script              //The script which runs GENIE3
    genie_con                //This script converts the GENIE3 output into something which can be measured
    num_genes               //This states the number of genes in the data


    //The nextflow process
    main:
    //STEP 1: run the Genie3 algorithm 
    genie_output = GENIE3_RUN(
        current,
        genie_script
    )
    //STEP2: Convert the file into a readbale format
    genie_matrix = GENIE_CONVERSION(
        genie_output,
        genie_con,
        num_genes
    )

    //We emit where the files are now located
    emit:
    genie_output = GENIE_CONVERSION.out
}

workflow{
    //Part 0: Using GENIE3 priors
    // genie_prior = GENIE3(
    //     messy,
    //     genie_script,
    //     genie_con, 
    //     num_genes
    // )


    //Part 1: Altering w0
    w0_TESTS(
        reads,
        original,
        p_val,
        input_fix,
        eb_script,
        prior_type, 
        keep,
        w0, 
        num_genes,
        original,
        data_type[0],
        method[0],
        num_cells,
        metrics,
        threshold_str,
        extra
    )

    //Part2: Altering prior formula
    PRIOR_CHANGE_TEST(
        reads,
        input_fix,
        prior_type,
        prior_eb,
        original,
        keep,
        multi,
        num_genes,
        original,
        data_type[0],
        method[1],
        num_cells,
        metrics,
        multi_str,
        p_val
    )
}


//GENIE3 prior testing:
//Runs GENIE
process GENIE3_RUN {

    cpus 4
    memory '4 GB'
    container 'genie3:latest'
    publishDir "${params.outdir}/genie3"

    input:
    path infile
    path genie_script

    output:
    path 'genie.txt'

    script:

    """
    Rscript ${genie_script} ${infile} > genie.txt
    """

} 

//Converts it into the workable form.
process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/genie3", mode: 'copy'

    input:
    path output_genie
    path genie_converter
    val num_genes

    output:
    path 'matrix_genie_normalised.csv'

    script:
    """
    python3 ${genie_converter} ${output_genie} ${num_genes}
    """
}