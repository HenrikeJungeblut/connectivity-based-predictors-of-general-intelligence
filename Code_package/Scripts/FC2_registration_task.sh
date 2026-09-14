# This script converts the Glasser (HCP-MMP1) surface atlas into a 3D volume 
# and warps it into native functional (fMRI) space. It then isolates the left 
# and right cortical hemispheres using intensity thresholds, and resamples them 
# to perfectly match the BOLD data's voxel grid using nearest-neighbor interpolation.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
  echo "registering $SUBJECT"
BASE_PATH="/path/to"
DERIV="$BASE_PATH/derivative/freesurfer/$SUBJECT/mri"
RAW="$BASE_PATH/raw_data/$SUBJECT/anat"

# Convert FreeSurfer annotation to T1 volume
mri_aparc2aseg \
  --s $SUBJECT \
  --annot HCP-MMP1 \
  --o $RAW/HCP-MMP1_T1.nii.gz

# apply transform created in FC1_registration_T1.sh
docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrtransform \
  "$RAW/HCP-MMP1_T1.nii.gz" \
  -linear "$RAW/${SUBJECT}_T1w2func_mrtrixformat.txt" \
  -inverse \
  "$DERIV/HCP-MMP1_task.nii.gz" \
  -force

# threshold to only retain cortex and split into right and left hemisphere
fslmaths $DERIV/HCP-MMP1_task.nii.gz \
  -thr 1001 -uthr 1180\
  $DERIV/HCP-MMP1_task_cortexonly_lh.nii.gz

fslmaths $DERIV/HCP-MMP1_task.nii.gz \
  -thr 2001 -uthr 2180\
  $DERIV/HCP-MMP1_task_cortexonly_rh.nii.gz

# resample to BOLD space
flirt -in $DERIV/HCP-MMP1_task_cortexonly_lh.nii.gz \
  -ref /path/to/derivative/$SUBJECT/tedana_switch/desc-denoised_bold.nii.gz \
  -applyxfm -usesqform \
  -out $DERIV/HCP-MMP1_task_cortexonly_lh_RS.nii.gz -interp nearestneighbour

flirt -in $DERIV/HCP-MMP1_task_cortexonly_rh.nii.gz \
  -ref /path/to/derivative/$SUBJECT/tedana_switch/desc-denoised_bold.nii.gz \
  -applyxfm -usesqform \
  -out $DERIV/HCP-MMP1_task_cortexonly_rh_RS.nii.gz -interp nearestneighbour

done

