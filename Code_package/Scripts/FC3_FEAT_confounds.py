# This script extracts 9 specific nuisance regressors (Global Signal, CSF, 
# White Matter, and 6 head motion parameters) from the fMRIPrep confounds 
# file. It replaces missing values with 0 and saves the data as a headerless 
# text file formatted specifically for use in FSL FEAT first-level analysis.

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd

for vp in ["sub-010"]:     # run for all subjects      
    # Load file
    confounds_file = pd.read_csv("/path/to/derivative/"+vp+ 
                                 "/fmriprep_switch/ses-01/func/"+vp+ 
                                 "_ses-01_task-switch_acq-mb3me4_desc-confounds_timeseries.tsv", sep='\t')

    # Extract head motion confounds, including temporal derivatives and quadratics
    motion_confounds = [x for x in confounds_file if 'trans' in x or 'rot' in x or 'global' in x or 'csf' in x or 'white' in x] 
    extracted_confounds = confounds_file[motion_confounds]
    extracted_confounds = extracted_confounds[['global_signal', 'csf', 'white_matter', 'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']]
    
    extracted_confounds = extracted_confounds.fillna(0)

    # Save confounds as .txt file
    extracted_confounds.to_csv("/path/to/derivative/" 
                               +vp+"/fmriprep_switch/ses-01/func/"+vp+ 
                               "_FEATconfounds.txt", sep='\t', header = False, index=False)
    
