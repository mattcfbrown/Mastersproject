//This script aims to take the threshold scripts and run them

//Scripts
file_correct = file( params.input_fix ) 
NI_script    = file( params.NI_script )
NI_con_scr   = file( params.NI_convert )
metric       = file( params.metrics )
genie_script = file( params.genie_script )
genie_con    = file( params.genie_con )
Emp_bay      = file( params.Emp_bay )


//Data
clean     = file( params.clean_data )
messy     = file( params.messy_data )
original  = file( params.original )
no_prior  = file( params.no_prior )

//Import the workflows needed
include { THRESH_INFORMATION_MEASURES as THRESH_INFORMATION_MEASURES_CLEAN } from './Modules'
include { THRESH_INFORMATION_MEASURES as THRESH_INFORMATION_MEASURES_MESSY } from './Modules'
include { THRESH_GENIE3 as THRESH_GENIE3_CLEAN } from './Modules'
include { THRESH_GENIE3 as THRESH_GENIE3_MESSY } from './Modules'
include { THRESH_BAYES as THRESH_BAYES_CLEAN_NO_PRIOR } from './Modules'
include { THRESH_BAYES as THRESH_BAYES_MESSY_NO_PRIOR } from './Modules'
include { THRESH_BAYES as THRESH_BAYES_CLEAN_PRIOR } from './Modules'
include { THRESH_BAYES as THRESH_BAYES_MESSY_PRIOR } from './Modules'


//Floating values
threshold_IM = [0.01, 0.05, 0.1, 0.2, 0.3]
threshold_str_IM = "0.01,0.05,0.1,0.2,0.3"
threshold_genie = [0.01,0.05,0.1,0.15,0.2]
threshold_str_genie = "0.01,0.05,0.1,0.15,0.2"
p_val = [0.5, 0.75, 0.9, 0.95, 0.99]
p_val_str = "0.5,0.75,0.9,0.95,0.99"
num_genes = 25
num_cells = 500
type = ["Clean", "Messy"]
method = ["InformationMeasures","GENIE3", "EB_NO_PRIOR", "EB_PRIOR"]
prior_type = ["zero", "full"]
keep = 0.9
gam_or_norm = "normal"
inference = "PUC"

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

    THRESH_BAYES_CLEAN_NO_PRIOR(
        clean,
        no_prior,
        p_val,
        file_correct,
        Emp_bay,
        prior_type[0],
        keep,
        gam_or_norm,
        inference,
        num_genes,
        original,
        type[0],
        method[2],
        num_cells,
        metric,
        p_val_str
    )

    THRESH_BAYES_MESSY_NO_PRIOR(
        messy,
        no_prior,
        p_val,
        file_correct,
        Emp_bay,
        prior_type[0],
        keep,
        gam_or_norm,
        inference,
        num_genes,
        original,
        type[1],
        method[2],
        num_cells,
        metric,
        p_val_str
    )

    THRESH_BAYES_CLEAN_PRIOR(
        clean,
        original,
        p_val,
        file_correct,
        Emp_bay,
        prior_type[1],
        keep,
        gam_or_norm,
        inference,
        num_genes,
        original,
        type[0],
        method[3],
        num_cells,
        metric,
        p_val_str
    )

    THRESH_BAYES_MESSY_PRIOR(
        messy,
        original,
        p_val,
        file_correct,
        Emp_bay,
        prior_type[1],
        keep,
        gam_or_norm,
        inference,
        num_genes,
        original,
        type[1],
        method[3],
        num_cells,
        metric,
        p_val_str
    )
}