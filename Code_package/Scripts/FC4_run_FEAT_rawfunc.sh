# This script performs nuisance regression to clean the fMRI data. It converts 
# the plain-text confound regressors into FSL's required matrix format (.mat) 
# and uses fsl_glm to regress them out of the denoised BOLD data, outputting 
# a clean, residualized 4D timeseries for connectivity analysis.

cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
# convert nuisance regressor txt to .mat
Text2Vest /path/to/derivative/$SUBJECT/fmriprep_switch/ses-01/func/${SUBJECT}_FEATconfounds.txt \
 /path/to/derivative/$SUBJECT/fmriprep_switch/ses-01/func/${SUBJECT}_FEATconfounds.mat

# run glm
fsl_glm -i /path/to/derivative/$SUBJECT/tedana_switch/desc-denoised_bold.nii.gz \
 -d /path/to/derivative/$SUBJECT/fmriprep_switch/ses-01/func/${SUBJECT}_FEATconfounds.mat \
  --out_res=/path/to/derivative/$SUBJECT/func/rawfunc_res4d.nii.gz
done

