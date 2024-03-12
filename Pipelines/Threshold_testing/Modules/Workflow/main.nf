//A script to hold all the individual workflows we wish to run



//Workflows
//Information measures
workflow INFORMATION_MEASURES {
    //Inputs
    take:
    current                     //Our original count matrix
    file_correct               //Puts our count matrix into a workable form
    NI_script                 //The script which runs Information measures
    NI_con_scr               //This script converts the Information measures output into something which can be measured
    threshold               //This is the threshold used
    num_genes              //The number of genes in the dataset


    //The nextflow process
    main:
    //STEP 1: Put the data into a readable format
    corrected_matrix = INPUT_INFORMATION_MEASURES(
        current,
        file_correct
    )
    //STEP 2: Perform the algorithm
    output_NI = INFORMATION_MEASURES_RUN(
        corrected_matrix,
        NI_script,
        threshold
    )
    //STEP 3: Convert it into a useable form
    ni_matrix = NI_CONVERSION(
        output_NI,
        NI_con_scr,
        num_genes,
        threshold
    )

    //We emit where the files are now located
    emit:
    im_output = NI_CONVERSION.out
}

//GENIE3 workflow
workflow GENIE3 {
    //Inputs
    take:
    current                    //Our original count matrix
    genie_script              //The script which runs GENIE3
    genie_con                //This script converts the GENIE3 output into something which can be measured
    threshold               //This is the threshold used
    num_genes              //This states the number of genes in the data


    //The nextflow process
    main:
    //STEP 1: run the Genie3 algorithm 
    genie_output = GENIE3_RUN_NOR(
        current,
        genie_script,
        threshold
    )
    //STEP2: Convert the file into a readbale format
    genie_matrix = GENIE_CONVERSION_NOR(
        genie_output,
        genie_con,
        num_genes,
        threshold
    )

    //We emit where the files are now located
    emit:
    genie_output = GENIE_CONVERSION_NOR.out
}

//Processes
//Mutual information process
//Part 1, convert into a useable form
process INPUT_INFORMATION_MEASURES {

    publishDir "${params.outdir}/Information_Measures"

    input:
    path infile
    path script

    output:
    path 'formatted_data.txt'

    script:

    """
    python3 ${script} ${infile} > formatted_data.txt
    """
}

//Part 2: Run the script
process INFORMATION_MEASURES_RUN {

    container 'networkinference:latest'
    publishDir "${params.outdir}/Information_Measures"

    input:
    path infile
    path jlscript
    val threshold

    output:
    path 'outfile_NI.txt'

    script:

    """

    julia ${jlscript} ${infile} ${threshold} > outfile_NI.txt

    """
}

//Part 3: Converted into something useful
process NI_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/Information_Measures"

    input:
    path output_NI
    path NI_converter
    val num_genes
    val threshold

    output:
    path "matrix_NI_${threshold}.csv"

    script:
    """
    python3 ${NI_converter} ${output_NI} ${num_genes} ${threshold}
    """
}

//GENIE3 process
//Part 1: Running GENIE
process GENIE3_RUN_NOR {

    cpus 4
    memory '4 GB'
    container 'genie3:latest'
    publishDir "${params.outdir}/genie3"

    input:
    path infile
    path genie_script
    val threshold

    output:
    path 'genie.txt'

    script:

    """
    Rscript ${genie_script} ${infile} ${threshold} > genie.txt
    """

}

//Part 2: Converting GENIE
process GENIE_CONVERSION_NOR {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/genie3", mode: 'copy'

    input:
    path output_genie
    path genie_converter
    val num_genes
    val threshold

    output:
    path "matrix_genie_${threshold}.csv"

    script:
    """
    python3 ${genie_converter} ${output_genie} ${num_genes} ${threshold}
    """
}
