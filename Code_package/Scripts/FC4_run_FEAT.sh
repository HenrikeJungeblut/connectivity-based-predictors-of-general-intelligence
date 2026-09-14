# This script automates FSL's FEAT for first-level fMRI analysis. It uses 'sed' 
# to replace placeholders in a template .fsf file with subject-specific paths 
# (denoised BOLD data, 3-column task timing EVs, and extracted confounds), and 
# then executes the FEAT analysis to generate the residualized timeseries.

#!/bin/bash
#this script runs the first level FEAT. 
FSLDATADIR=/path/to/derivative
COLUMNFILESDIR=/path/to/derivative/scripts/2_scripts_BOLD/3ColumnFiles
 
#########
 
#makes the orient file
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
for sub in 021; do
 mkdir ${FSLDATADIR}/sub-${sub}/func/fsf_file_residualization
 OUTPUT=${FSLDATADIR}/sub-${sub}/func/first_level_residual.feat 
 DATA=${FSLDATADIR}/sub-${sub}/tedana_switch/desc-denoised_bold.nii.gz
 EV1=${COLUMNFILESDIR}/S0${sub}_LM.txt
 EV2=${COLUMNFILESDIR}/S0${sub}_OE.txt
 EV3=${COLUMNFILESDIR}/S0${sub}_SW.txt
 CONFOUND=${FSLDATADIR}/sub-${sub}/fmriprep_switch/ses-01/func/sub-${sub}_FEATconfounds.txt
 echo $OUTPUT
 #makes the fsf files from the template fsf file
 for i in 'template.fsf'; do
  sed -e 's@OUTPUT@'$OUTPUT'@g' \
   -e 's@ANAT@'$ANAT'@g' \
   -e 's@EV1@'$EV1'@g' \
   -e 's@EV2@'$EV2'@g' \
   -e 's@EV3@'$EV3'@g' \
   -e 's@CONFOUND@'$CONFOUND'@g' \
   -e 's@DATA@'$DATA'@g' <$i> ${FSLDATADIR}/sub-${sub}/func/fsf_file_residualization/FEAT_sub-${sub}.fsf 
 done
 #runs the analysis using the newly created fsf file
 feat ${FSLDATADIR}/sub-${sub}/func/fsf_file_residualization/FEAT_sub-${sub}.fsf
done
