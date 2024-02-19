//All the scripts required are listed below
nlnet_rscript  = file( params.nlnet_rscript )
nlnet_con_scr  = file( params.nlnet_convert )
file_correct   = file( params.input_fix )
NI_script      = file( params.NI_script )
NI_con_scr     = file( params.NI_convert )
genie_script   = file( params.genie3_rscript )
genie_con      = file( params.genie_convert )
slingshot_scr  = file( params.slingshot )
SCODE_scr      = file( params.SCODE )
scode_con_scr  = file( params.SCODE_convert )
metric_program = file( params.metric )

//All the data which will be used is listed below
gene_read10    = file( params.gene_read10 )
gene_count10   = file( params.gene_count10 )
gene_orig10    = file( params.gene_orig10 )
gene_read25    = file( params.gene_read25 )
gene_count25   = file( params.gene_count25 )
gene_orig25    = file( params.gene_orig25 )

//Now get the packages we wish to work with
include { MAJOR as TEN_GENES } from './modules'
include { MAJOR as TWENTY_FIVE_GENES } from './modules'

threshold = 0.1

workflow {

    //Firstly run for 10 genes
    TEN_GENES(
        gene_read10,
        nlnet_rscript,
        nlnet_con_scr,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold,
        genie_script,
        genie_con,
        slingshot_scr,
        scode_con_scr,
        gene_count10,
        gene_orig10,
        metric_program,
        10,
        SCODE_scr,
        100
    )

    //25 genes
    TWENTY_FIVE_GENES(
        gene_read25,
        nlnet_rscript,
        nlnet_con_scr,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold,
        genie_script,
        genie_con,
        slingshot_scr,
        scode_con_scr,
        gene_count25,
        gene_orig25,
        metric_program,
        25,
        SCODE_scr,
        2000
    )
}