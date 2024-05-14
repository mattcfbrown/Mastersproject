// This here will aim to do the following:
    //It will take a section of the read data and spit out a belief on how the network is formed
    //It will do this several times, depending on the size of the network, from here you will have several guesses
    //Can make a consensus, you decide what the network should be based on what you see
//

//Scripts used
convert_file = file( params.convert_file )

//Data used
clean = file( params.clean )

//Floating variables created
num_cells_taken = 2000

//Workflows used


//Workflow where we run stuff
workflow {
    //Here is the data which will be used
    ch_data = Channel.of(
        tuple(clean,convert_file,1,num_cells_taken),
        tuple(clean,convert_file,2,num_cells_taken),
        tuple(clean,convert_file,3,num_cells_taken),
        tuple(clean,convert_file,4,num_cells_taken)
    )

    ch_data | EB_CONVERT

}

//Processes used
//This converts the EmpiricalBayes script into a useable form.
process EB_CONVERT{

    container 'nlnet_convert:latest'
    publishDir "${params.outdir}/Converted_files"

    input:
    tuple path(infile), path(script), val(iter_no), val(num_cells_taken)

    output:
    path "formatted_data_${iter_no}.txt"

    script:

    """
    python3 ${script} ${infile} ${num_cells_taken} ${iter_no}
    """    

}