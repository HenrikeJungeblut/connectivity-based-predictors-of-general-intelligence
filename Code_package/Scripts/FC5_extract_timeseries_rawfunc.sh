# This script extracts the mean BOLD time series for each cortical region of 
# the Glasser atlas from the cleaned fMRI (output of FC4_run_FEAT_rawfunc) data using FSL's fslmeants. It 
# processes the left and right hemispheres separately, and then uses awk to 
# filter out any zero-valued columns (empty/missing regions) from the outputs.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
# extract timeseries
fslmeants -i /path/to/derivative/$SUBJECT/func/rawfunc_res4d.nii.gz \
 -o /path/to/derivative/$SUBJECT/tedana_switch/rawfunc_${SUBJECT}_lh_timeseries_per_glasser_region.txt \
 --label=/path/to/derivative/freesurfer/$SUBJECT/mri/HCP-MMP1_task_cortexonly_lh_RS.nii.gz

fslmeants -i /path/to/derivative/$SUBJECT/func/rawfunc_res4d.nii.gz \
 -o /path/to/derivative/$SUBJECT/tedana_switch/rawfunc_${SUBJECT}_rh_timeseries_per_glasser_region.txt \
 --label=/path/to/derivative/freesurfer/$SUBJECT/mri/HCP-MMP1_task_cortexonly_rh_RS.nii.gz

#remove columns with only zeros
awk '{for(i=1; i<=NF; i++) if($i != 0) printf "%s ", $i; print ""}' /path/to/derivative/$SUBJECT/tedana_switch/rawfunc_${SUBJECT}_lh_timeseries_per_glasser_region.txt > /path/to/MRI/func/rawfunc_timeseries_per_region/rawfunc_${SUBJECT}_lh_filtered_timeseries_per_glasser_region.txt
awk '{for(i=1; i<=NF; i++) if($i != 0) printf "%s ", $i; print ""}' /path/to/derivative/$SUBJECT/tedana_switch/rawfunc_${SUBJECT}_rh_timeseries_per_glasser_region.txt > /path/to/MRI/func/rawfunc_timeseries_per_region/rawfunc_${SUBJECT}_rh_filtered_timeseries_per_glasser_region.txt
done

