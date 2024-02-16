

// declare nlnet rscript IS a file (so we can pass it correctly to process)
nlnet_rscript  = file( params.nlnet_rscript )
NI_script      = file( params.NI_script )
EB_script      = file( params.EB_script )
count_matrix   = file( params.count_data )
count_EB       = file( params.EB_count )
SCENIC_test    = file( params.test_SC)
SCENIC_matrix  = file( params.matrix)
Trans_fact     = file( params.TFs)
Network        = file( params.network)
toy            = file( params.toy)
toy_original   = file( params.toy_original)
nlnet_con_scr  = file( params.nlnet_convert)
original_mat   = file( params.original_mat)
metric_program = file( params.metric)
NI_con_scr     = file( params.NI_convert)
file_correct   = file( params.input_fix)
genie_script   = file( params.genie3_rscript)
genie_con      = file( params.genie_convert)
slingshot_scr  = file( params.slingshot)
toy_counts     = file( params.toy_count)
SCODE_scr      = file( params.SCODE)
scode_con_scr  = file( params.SCODE_convert)


//Data we are currently using
current   = toy
original  = toy_original
threshold = 0.1
num_genes = 9
num_cells = 300

//LNET workflow
workflow NLNET {
    //Inputs
    take:
    path current                    //Our original count matrix
    path nlnet_rscript             //The script which runs NLNET
    path nlnet_con_scr            //This script converts the NLNET output into something which can be measured
    val num_genes                //This is the number of genes being used

    //The nextflow process
    main:
    //STEP 1: Run NLNET
    output_nlnet = NLNET_RUN(
        current,
        nlnet_rscript
    )
    //STEP 2: Convert into a useable form 
    nlnet_matrix = NLNET_CONVERSION (
        output_nlnet,
        nlnet_con_scr,
        num_genes
    )

    //We emit where the files are now located
    emit:
    nlnet_output = NLNET_CONVERSION.out
}

//Information measures workflow
workflow INFORMATION_MEASURES {
    //Inputs
    take:
    path current                    //Our original count matrix
    path file_correct              //Puts our count matrix into a workable form
    path NI_script                //The script which runs Information measures
    path NI_con_scr              //This script converts the Information measures output into something which can be measured
    val threshold               //This is the threshold used


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
        num_genes
    )

    //We emit where the files are now located
    emit:
    im_output = NI_CONVERSION.out
}

//GENIE3 Workflow
workflow GENIE3 {
    //Inputs
    take:
    path current                   //Our original count matrix
    path genie_script             //The script which runs GENIE3
    path genie_con               //This script converts the GENIE3 output into something which can be measured
    val threshold               //This is the threshold used
    val num_genes              //This states the number of genes in the data


    //The nextflow process
    main:
    //STEP 1: run the Genie3 algorithm 
    genie_output = GENIE3_RUN(
        current,
        genie_script,
        threshold
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

//SCODE Workflow
workflow SCODE{
    //Inputs
    take:
    path slingshot_scr                  //The script which runs the slingshot algorithm
    path toy_counts                    //A matrix displaying the counts of each gene in the individual experiment
    path current                      //Our original count matrix
    path SCODE_scr                   //The script which runs SCODE
    path scode_con_scr              //This script converts the SCODE output into something which can be measured
    val num_cells                  //This is the number of cells in the data
    val num_genes                 //This states the number of genes in the data    

    main:
    //Step 1: Get pseudotime inputs using slingshot
    time = SLINGSHOT(
        slingshot_scr,
        toy_counts
    )
    //Step 2: Run the SCODE code
    output_SCODE = SCODE_RUN(
        SCODE_scr,
        current,
        time,
        num_genes,
        num_cells       
    )
    //Step 3: Convert the output into something workable
    SCODE_CONVERT(
        scode_con_scr,
        output_SCODE
    )   

    //Output of the workflow
    emit:
    SCODE_output = SCODE_CONVERT.out
}

workflow {
    //NLNET workflow
    NLNET(
        current,
        nlnet_rscript,
    )

    //Information measures workflow
    INFORMATION_MEASURES(
        current,
        file_correct
    )

    //GENIE3 Method
    GENIE3(
        current,
        genie_script
    )

    //SCODE Method
    SCODE(
        slingshot_scr,
        toy_counts        
    )

    //Metrics workflow
    METRICS(
        NLNET.out.nlnet_output,
        INFORMATION_MEASURES.out.im_output,
        GENIE3.out.genie_output,
        SCODE.out.SCODE_output,
        original,
        metric_program
    )

    //Empirical Bayes method for network inference
    //EMPIRICAL_BAYES(
    //    count_EB,
    //    EB_script
    //).view()

    //SCENIC(
    //    SCENIC_matrix,
    //    Trans_fact
    //).view()


}

process NLNET_RUN {

    container 'nlnet:latest'
    publishDir "${params.outdir}/nlnet"

    input:
    path infile
    path rscript

    output:
    path 'outfile_nlnet.txt'

    script:
    """
    Rscript ${rscript} ${infile} > outfile_nlnet.txt
    """
}

process NLNET_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/nlnet"

    input:
    path output_nlnet
    path nlnet_converter
    val num_genes

    output:
    path 'matrix_nlnet.csv'

    script:
    """
    python3 ${nlnet_converter} ${output_nlnet} ${num_genes} > matrix_nlnet.csv
    """

}

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

process NI_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/Information_Measures"

    input:
    path output_NI
    path NI_converter
    val num_genes

    output:
    path 'matrix_NI.csv'

    script:
    """
    python3 ${NI_converter} ${output_NI} ${num_genes} > matrix_NI.csv
    """
}

process GENIE3_RUN {

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

process GENIE_CONVERSION {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/genie3"

    input:
    path output_genie
    path genie_converter
    val num_genes

    output:
    path 'matrix_genie.csv'

    script:
    """
    python3 ${genie_converter} ${output_genie} ${num_genes} > matrix_genie.csv
    """
}

process SLINGSHOT {

    container 'slingshot:latest'
    publishDir "${params.outdir}/SCODE"

    input:
    path slingshot_script
    path toy_counts

    output:
    path 'mytable_R.txt'

    script:
    """
    Rscript ${slingshot_script} ${toy_counts} > mytable_R.txt
    """
}

process SCODE_RUN {

    publishDir "${params.outdir}/SCODE"

    input:
    path SCODE_script
    path expression_data
    path pseudotime
    val num_genes
    val num_cells

    output:
    path './out' 

    script:
    """
    Rscript ${SCODE_script} ${expression_data} ${pseudotime} out ${num_genes} 4 ${num_cells} 100
    """
}

process SCODE_CONVERT {

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/SCODE"

    input:
    path scode_con_scr
    path input
    
    output:
    path 'matrix_SCODE.csv'

    script:
    """
    python3 ${scode_con_scr} ${input}/A.txt > matrix_SCODE.csv
    """
}

process EMPIRICAL_BAYES {
    
    container 'empiricalbayes:latest'

    input:
    path infile
    path jlscript

    output:
    path 'outfile_EB.txt'

    script:

    """

    julia ${jlscript} ${infile} > outfile_EB.txt

    """
}

process SCENIC {

    container 'aertslab/pyscenic:0.9.18'
    publishDir "${params.outdir}/SCENIC"
    
    input:
    path grn
    path TF
    
    output:
    path "3A_seurat_RNA_pySCENIC_GRN_adjacencies.csv", emit: outfile
    
    script:

    """
    
    pyscenic grn ${grn} ${TF} -o expr_mat.adjacencies.tsv


    """
}

process METRICS {

    container 'metrics:latest'
    publishDir "${params.outdir}/metrics"

    input:
    path nlnet_matrix
    path ni_matrix
    path genie_matrix
    path scode_matrix
    path original_matrix
    path metric_script

    output:
    path 'ROCplot.pdf'
    path 'matthew.txt'

    script:

    """
    python3 ${metric_script} ${nlnet_matrix} ${ni_matrix} ${genie_matrix} ${scode_matrix} ${original_matrix}
    """    

}