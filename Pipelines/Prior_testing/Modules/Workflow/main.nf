//This here will contain all the functions which we will be using


//Empirical Bayes workflow
workflow EMPIRICAL_BAYES{
    take:
    reads                              //This is the gene sequencing reads
    priors                            //This is the prior information 
    p_val                            //This is the p value
    conversion_script               //Script which converts the files into something workable
    eb_script                      //This script runs Empirical Bayes
    type                          //This is the type of prior we are giving it
    keep                         //The proportion of values we are keeping when running EB
    w0                          //The parameter used to calculate the prior information
    extra                      // So I don't have to edit my code

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
        w0,
        extra
    )

    emit:
    output = EMPIRICAL_BAYES_RUN.out
}

//Prior change workflow
workflow PRIOR_CHANGE{
    take:
    reads                                           //The data we are analysing
    conversion_script                              //Converts the data into something useable
    prior_type                                    //The type of prior information we are using
    prior_script                                 //The script which runs out EB test
    priors                                      //The prior information used
    to_keep                                    //The percentage of the values to keep
    multi                                     //How we are affecting the priors
    p_val                                    //P value used

    main:
    //Step 1: Put it in a workable form
    data = EB_CONVERT(
        reads,
        conversion_script,
        prior_type
    )

    //Step 2: Now run the prior script
    EMPIRICAL_BAYES_PRIOR(
        prior_script,
        data,
        priors,
        to_keep,
        multi,
        p_val,
        prior_type
    )

    emit:
    EMPIRICAL_BAYES_PRIOR.out
}

//Threshold Workflow
workflow THRESHOLD{

    take:
    metric
    threshold1
    threshold2
    threshold3    
    threshold4
    threshold5
    original
    threshold_str
    type
    method
    num_cells
    num_genes

    main:
    METRIC_THRESHOLD(
        metric,
        threshold1,
        threshold2,
        threshold3,   
        threshold4,
        threshold5,
        original,
        threshold_str,
        type,
        method,
        num_cells,
        num_genes
    )

}

//Empirical Bayes process
//Step 1: Convert it into a useable form
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

//This gets the output of the Empirical bayes run
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
    val w0
    val extra

    output:
    path "Eb_matrix_${type}_${w0}_${extra}.csv"

    script:

    """
    julia ${script} ${data} ${p_val} ${priors} ${type} ${keep} ${w0} ${extra}
    """
}

//Prior formula:
process EMPIRICAL_BAYES_PRIOR{

    container 'changing_priors:latest'
    publishDir "${params.outdir}/multi"

    input:
    path script
    path data
    path priors
    val to_keep
    val multi
    val p_val
    val type


    output:
    path "Eb_matrix_${multi}_${type}.csv"

    script:

    """
    julia ${script} ${data} ${priors} ${to_keep} ${multi} ${p_val} ${type}
    """    
}

//Processes
//Metrics
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