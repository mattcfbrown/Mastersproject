//This script aims to take the threshold scripts and run them

//Scripts
file_correct = file( params.input_fix ) 
NI_script    = file( params.NI_script )
NI_con_scr   = file( params.NI_convert )
metric       = file( params.metrics )
genie_script = file( params.genie_script )
genie_con    = file( params.genie_con )

//Data
clean     = file( params.clean_data )
messy     = file( params.messy_data )
original  = file( params.original )

//Import the workflows needed
include { THRESH_INFORMATION_MEASURES as THRESH_INFORMATION_MEASURES_CLEAN } from './Modules'
include { THRESH_INFORMATION_MEASURES as THRESH_INFORMATION_MEASURES_MESSY } from './Modules'
include { THRESH_GENIE3 as THRESH_GENIE3_CLEAN } from './Modules'
include { THRESH_GENIE3 as THRESH_GENIE3_MESSY } from './Modules'

//Floating values
threshold_IM = [0.01, 0.05, 0.1, 0.2, 0.3]
threshold_str_IM = "0.01,0.05,0.1,0.2,0.3"
threshold_genie = [0.01,0.05,0.1,0.15,0.2]
threshold_str_genie = "0.01,0.05,0.1,0.15,0.2"
num_genes = 25
num_cells = 500
type = ["Clean", "Messy"]
method = ["InformationMeasures","GENIE3"]

//Main workflow
workflow{

    THRESH_INFORMATION_MEASURES_CLEAN(
        clean,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold_IM,
        num_genes,
        original,
        type[0],
        method[0],
        num_cells,
        metric,
        threshold_str_IM
    )

    THRESH_INFORMATION_MEASURES_MESSY(
        messy,
        file_correct,
        NI_script,
        NI_con_scr,
        threshold_IM,
        num_genes,
        original,
        type[1],
        method[0],
        num_cells,
        metric,
        threshold_str_IM        
    )

    THRESH_GENIE3_CLEAN(
        clean,
        genie_script,
        genie_con,
        threshold_genie,
        num_genes,
        original,
        type[0],
        method[1],
        num_cells,
        metric,
        threshold_str_genie
    )

    THRESH_GENIE3_MESSY(
        messy,
        genie_script,
        genie_con,
        threshold_genie,
        num_genes,
        original,
        type[1],
        method[1],
        num_cells,
        metric,
        threshold_str_genie
    )
}