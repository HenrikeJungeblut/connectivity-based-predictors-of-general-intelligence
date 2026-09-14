# This script uses FSL's eddy to correct the raw DWI data for eddy current-
# induced distortions, subject head motion, and susceptibility distortions (by 
# applying the topup results). It also uses the --repol flag to automatically 
# detect and replace slice-level signal dropouts (outliers).
cat /path/to/derivative/subject_names.txt | tr -d '\r' | while read sub; do

eddy --imain=/path/to/raw_data/${sub}/dwi/${sub}_dwi_AP.nii.gz \
--mask=/path/to/derivative/${sub}/dwi/hifi_nodif_brain_mask.nii.gz \
--bvals=/path/to/derivative/${sub}/dwi/bvals \
--bvecs=/path/to/derivative/${sub}/dwi/bvecs \
--acqp=/path/to/derivative/2_scripts_DWI/acq_param.txt \
--index=/path/to/derivative/2_scripts_DWI/index.txt \
--out=/path/to/derivative/${sub}/dwi/eddy_corrected_data \
--topup=/path/to/derivative/${sub}/dwi/topup_AP_PA_b0 --repol --data_is_shelled \
--verbose

echo 'Done SUB-'$sub
done


