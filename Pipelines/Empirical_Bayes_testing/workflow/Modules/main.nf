//This here is a program which is able to test how effective prior usage is for different softwares


//All the scripts needed
file_correct = file( params.input_fix )
bayes_script = file( params.Emp_Bayes )

//All the datasets needed
test_data    = file( params.test_data )
gene_read25  = file( params.gene_read25 )
gene_orig25  = file( params.gene_orig25)

//All floating variables
p_value = 0.9

//Workflow to run the whole thing

workflow EMPIRICAL_BAYES{
    take:
    reads                              //This is the gene sequencing reads
    priors                            //This is the prior information 
    p_val                            //This is the p value
    conversion_script               //Script which converts the files into something workable
    eb_script                      //This script runs Empirical Bayes
    type                          //This is the type of prior we are giving it
    keep                         //The proportion of values we are keeping when running EB
    gam_or_norm                 //This states whether to use Gamma or Normal fitting
    inference                  //What type of inference to use

    main:
    //Step 1: Makes it readable
    data = EB_CONVERT(
        reads,
        conversion_script,
        type
    )
    //Step 2: Run Empirical Bayes
    EMPIRICAL_BAYES_RUN(
        data,
        eb_script,
        priors,
        p_val,
        type,
        keep,
        gam_or_norm,
        inference
    )

    emit:
    output = EMPIRICAL_BAYES_RUN.out
}

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

//Workflow
workflow {

    //EMPIRICAL_BAYES(
    //    gene_read25,
    //    gene_orig25,
    //    p_value,
    //    file_correct,
    //    bayes_script
    //).view()


}

//This process converts the reads into something empirical Bayes can use
process EB_CONVERT{

    publishDir "${params.outdir}/${type}"

    input:
    path infile
    path script
    val type

    output:
    path 'formatted_data.txt'

    script:

    """
    python3 ${script} ${infile} > formatted_data.txt
    """    

}

process EMPIRICAL_BAYES_RUN {

    container 'empiricalbayes:latest'
    publishDir "${params.outdir}/${type}"

    input:
    path data
    path script
    path priors
    val p_val
    val type
    val keep
    val gam_or_norm
    val inference

    output:
    path "Eb_matrix_${type}.csv"

    script:

    """
    julia ${script} ${data} ${p_val} ${priors} ${type} ${keep} ${gam_or_norm} ${inference}
    """
   
}

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

process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/genie3", mode: 'copy'

    input:
    path output_genie
    path genie_converter
    val num_genes

    output:
    path 'matrix_genie.csv'

    script:
    """
    python3 ${genie_converter} ${output_genie} ${num_genes}
    """
}
